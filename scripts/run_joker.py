# Standard library
from os import path
import os
import time

# Third-party
import astropy.units as u
import h5py
import numpy as np
from schwimmbad import choose_pool
from thejoker.log import log as joker_logger
from thejoker.sampler import TheJoker, JokerParams

# Project
from twoface import unimodal_P
from twoface.log import log as logger
from twoface.sample_prior import make_prior_cache

from apogeebh.db import db_connect
from apogeebh.db import AllStar, StarResult, Status


def main(pool, seed, overwrite=False, _continue=False):

    # HACK: hard-coded configuration!
    db_path = path.abspath('../cache/apogeebh.sqlite')
    prior_samples_file = path.abspath('../cache/prior-samples.hdf5')
    results_filename = path.abspath('../cache/apogeebh-joker.hdf5')
    n_prior = 536870912 # number of prior samples to generate
    n_requested_samples = 256 # how many samples to generate, nominally
    max_samples_per_star = 2048 # max. number of posterior samples to save
    P_min = 1 * u.day
    P_max = 1024 * u.day
    jitter = 150 * u.m/u.s

    if not os.path.exists(db_path):
        raise IOError("sqlite database not found at '{0}'\n Did you run "
                      "scripts/load_dr15_db.py yet for that database?"
                      .format(db_path))

    logger.debug("Connecting to sqlite database at '{0}'".format(db_path))
    Session, engine = db_connect(database_path=db_path,
                                 ensure_db_exists=False)
    session = Session()

    # Retrieve or create a JokerRun instance
    params = JokerParams(P_min=P_min, P_max=P_max, jitter=jitter)

    # Create TheJoker sampler instance with the specified random seed and pool
    rnd = np.random.RandomState(seed=seed)
    logger.debug("Creating TheJoker instance with {0}, {1}".format(rnd, pool))
    joker = TheJoker(params, random_state=rnd, pool=pool)

    # Create a cache of prior samples (if it doesn't exist) and store the
    # filename in the database.
    if not os.path.exists(prior_samples_file) or overwrite:
        logger.debug("Prior samples file not found - generating {0} samples..."
                     .format(n_prior))
        make_prior_cache(prior_samples_file, joker,
                         nsamples=n_prior)
        logger.debug("...done")

    # Get done APOGEE ID's
    done_subq = session.query(AllStar.apogee_id)\
                       .join(StarResult, Status)\
                       .filter(Status.id > 0).distinct()

    # Query to get all stars associated with this run that need processing:
    # they should have a status id = 0 (needs processing)
    star_query = session.query(AllStar)\
                        .filter(AllStar.vscatter >= 5.)\
                        .filter(~AllStar.apogee_id.in_(done_subq))

    # Base query to get a StarResult for a given Star so we can update the
    # status, etc.
    result_query = session.query(StarResult).join(AllStar)\
                          .filter(Status.id == 0)\
                          .filter(~AllStar.apogee_id.in_(done_subq))

    n_stars = star_query.count()
    logger.info("{0} stars left to process".format(n_stars))

    # Ensure that the results file exists - this is where we cache samples that
    # pass the rejection sampling step
    if not os.path.exists(results_filename):
        with h5py.File(results_filename, 'w') as f:
            pass

    # --------------------------------------------------------------------------
    # Here is where we do the actual processing of the data for each star. We
    # loop through all stars that still need processing and iteratively
    # rejection sample with larger and larger prior sample batch sizes. We do
    # this for efficiency, but the argument for this is somewhat made up...

    count = 0 # how many stars we've processed in this star batch
    batch_size = 16 # MAGIC NUMBER: how many stars to process before committing
    for star in star_query.all():

        if result_query.filter(AllStar.apogee_id == star.apogee_id).count() < 1:
            logger.debug('Star {0} has no result object!'
                         .format(star.apogee_id))
            result = StarResult()
            star.result = result
            session.add(result)
            session.commit()

        # Retrieve existing StarResult from database. We limit(1) because
        # the APOGEE_ID isn't unique, but we attach all visits for a given
        # star to all rows, so grabbing one of them is fine.
        result = result_query.filter(AllStar.apogee_id == star.apogee_id)\
                             .limit(1).one()

        logger.log(1, "Starting star {0}".format(star.apogee_id))
        logger.log(1, "Current status: {0}".format(str(result.status)))
        t0 = time.time()

        data = star.apogeervdata()
        logger.log(1, "\t visits loaded ({:.2f} seconds)"
                   .format(time.time()-t0))
        try:
            samples, ln_prior = joker.iterative_rejection_sample(
                data=data, n_requested_samples=n_requested_samples,
                prior_cache_file=prior_samples_file,
                n_prior_samples=n_prior, return_logprobs=True)

        except Exception as e:
            logger.warning("\t Failed sampling for star {0} \n Error: {1}"
                           .format(star.apogee_id, str(e)))
            continue

        logger.debug("\t done sampling ({:.2f} seconds)".format(time.time()-t0))

        # For now, it's sufficient to write the run results to an HDF5 file
        all_ln_probs = ln_prior[:max_samples_per_star]
        samples = samples[:max_samples_per_star]
        n_actual_samples = len(all_ln_probs)

        # Write the samples that pass to the results file
        with h5py.File(results_filename, 'r+') as f:
            if star.apogee_id in f:
                del f[star.apogee_id]

            # HACK: this will overwrite the past samples!
            g = f.create_group(star.apogee_id)
            samples.to_hdf5(g)

            if 'ln_prior_probs' in g:
                del g['ln_prior_probs']
            g.create_dataset('ln_prior_probs', data=all_ln_probs)

        logger.debug("\t saved samples ({:.2f} seconds)".format(time.time()-t0))

        if n_actual_samples >= n_requested_samples:
            result.status_id = 4 # completed

        elif n_actual_samples == 1:
            # Only one sample was returned - this is probably unimodal, so this
            # star needs MCMC
            result.status_id = 2 # needs mcmc

        else:

            if unimodal_P(samples, data):
                # Multiple samples were returned, but they look unimodal
                result.status_id = 2 # needs mcmc

            else:
                # Multiple samples were returned, but not enough to satisfy the
                # number requested in the config file
                result.status_id = 1 # needs more samples

        logger.debug("...done with star {} ({:.2f} seconds)"
                     .format(star.apogee_id, time.time()-t0))

        if count % batch_size == 0 and count > 0:
            session.commit()

        count += 1

    pool.close()

    session.commit()
    session.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")

    vq_group = parser.add_mutually_exclusive_group()
    vq_group.add_argument('-v', '--verbose', action='count', default=0,
                          dest='verbosity')
    vq_group.add_argument('-q', '--quiet', action='count', default=0,
                          dest='quietness')

    parser.add_argument("--overwrite", dest="overwrite", default=False,
                        action="store_true",
                        help="Overwrite any existing results.")

    parser.add_argument("-s", "--seed", dest="seed", default=None, type=int,
                        help="Random number seed")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--procs", dest="n_procs", default=1,
                       type=int, help="Number of processes.")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    args = parser.parse_args()

    loggers = [joker_logger, logger]

    # Set logger level based on verbose flags
    if args.verbosity != 0:
        if args.verbosity == 1:
            logger.setLevel(logging.DEBUG)
        else: # anything >= 2
            logger.setLevel(1)
            joker_logger.setLevel(1)

    elif args.quietness != 0:
        if args.quietness == 1:
            logger.setLevel(logging.WARNING)
            joker_logger.setLevel(logging.WARNING)
        else: # anything >= 2
            logger.setLevel(logging.ERROR)
            joker_logger.setLevel(logging.ERROR)

    else: # default
        logger.setLevel(logging.INFO)
        joker_logger.setLevel(logging.INFO)

    pool_kwargs = dict(mpi=args.mpi, processes=args.n_procs)
    pool = choose_pool(**pool_kwargs)

    main(pool=pool, seed=args.seed, overwrite=args.overwrite)
