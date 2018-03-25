"""
Generate posterior samples over orbital parameters for all data files in a
specified directory.
"""

# Standard library
import glob
from os import path
import os
import pickle
import sys
import time

# Third-party
from astropy.table import QTable
import astropy.units as u
import h5py
import numpy as np
from schwimmbad import choose_pool
from thejoker.log import log as joker_logger
from thejoker import TheJoker, RVData, JokerSamples
from thejoker.sampler.mcmc import TheJokerMCMCModel
import yaml
import emcee

# Project
from twoface.log import log as logger
from twoface.sample_prior import make_prior_cache
from twoface.util import config_to_jokerparams
from twoface.samples_analysis import unimodal_P, MAP_sample
from twoface.mcmc_helpers import gelman_rubin


global data, joker_params


def logprob(p):
    model = TheJokerMCMCModel(joker_params, data)
    return model(p)


def main(data_path, apogee_id, config, data_file_ext, pool, overwrite=False):

    cache_path = path.join(data_path, 'cache')
    os.makedirs(cache_path, exist_ok=True)

    # Create TheJoker sampler instance with the specified random seed and pool
    rnd = np.random.RandomState()
    logger.debug("Creating TheJoker instance with {0}, {1}".format(rnd, pool))
    joker = TheJoker(joker_params, random_state=rnd, pool=pool)

    n_walkers = config['emcee']['n_walkers']
    # n_steps = config['emcee']['n_steps']
    n_steps = 32

    mcmc_model_filename = path.join(cache_path, 'model.pickle')

    logger.info('Processing {0}'.format(apogee_id))

    joker_results_filename = path.join(cache_path,
                                       '{0}-joker.hdf5'.format(apogee_id))
    mcmc_results_filename = path.join(cache_path,
                                      '{0}-mcmc.hdf5'.format(apogee_id))
    mcmc_chain_filename = path.join(cache_path,
                                    '{0}-chain.npy'.format(apogee_id))

    with h5py.File(joker_results_filename) as f:
        joker_samples = JokerSamples.from_hdf5(f)

    if not path.exists(mcmc_chain_filename) or overwrite:
        joker = TheJoker(joker_params)
        joker.params.jitter = (8.5, 0.9) # HACK!

        model = TheJokerMCMCModel(joker_params=joker_params, data=data)

        if unimodal_P(joker_samples, data):
            logger.debug("Samples are unimodal. Preparing to run MCMC...")

            sample = MAP_sample(data, joker_samples, joker.params)

            ball_scale = 1E-5
            p0_mean = np.squeeze(model.pack_samples(sample))

            # P, M0, e, omega, jitter, K, v0
            p0 = np.zeros((n_walkers, len(p0_mean)))
            for i in range(p0.shape[1]):
                if i in [2, 4]: # eccentricity, jitter
                    p0[:, i] = np.abs(np.random.normal(p0_mean[i], ball_scale,
                                                       size=n_walkers))

                else:
                    p0[:, i] = np.random.normal(p0_mean[i], ball_scale,
                                                size=n_walkers)

            p0 = model.to_mcmc_params(p0.T).T

            n_dim = p0.shape[1]
            sampler = emcee.EnsembleSampler(n_walkers, n_dim, model,
                                            pool=pool)

            logger.debug('Running MCMC for {0} steps...'.format(n_steps))
            time0 = time.time()
            _ = sampler.run_mcmc(p0, n_steps)
            logger.debug('...time spent sampling: {0}'
                         .format(time.time() - time0))

            samples = model.unpack_samples_mcmc(sampler.chain[:, -1])
            samples.t0 = joker_samples.t0

            np.save(mcmc_chain_filename, sampler.chain.astype('f4'))

            if not path.exists(mcmc_model_filename):
                with open(mcmc_model_filename, 'wb') as f:
                    pickle.dump(model, f)

        else:
            logger.error("Samples are multimodal!")
            return

    if not path.exists(mcmc_results_filename) or overwrite:
        with h5py.File(mcmc_results_filename) as f:
            chain = np.load(mcmc_chain_filename)
            n_walkers, n_steps, n_pars = chain.shape

            g = f.create_group(apogee_id)

            logger.debug('Adding star {0} to MCMC cache'.format(star.apogee_id))
            try:
                g2 = g.create_group('chain-stats')
                Rs = gelman_rubin(chain[:, n_steps//2:])
                g2.create_dataset(name='gelman_rubin', data=Rs)

                # take the last sample, downsample
                end_pos = chain[:run.requested_samples_per_star, -1]
                samples = model.unpack_samples_mcmc(end_pos)
                samples.to_hdf5(g)

            except Exception as e:
                raise

            finally:
                del g


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

    parser.add_argument("-o", "--overwrite", dest="overwrite", default=False,
                        action='store_true')

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--procs", dest="n_procs", default=1,
                       type=int, help="Number of processes.")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    parser.add_argument("--apogeeid", dest="apogee_id", required=True,
                        type=str, help="The source to run on.")
    parser.add_argument("--data-path", dest="data_path", required=True,
                        type=str, help="Path to the data files to run on.")
    parser.add_argument("-c", "--config", dest="config_file", required=True,
                        type=str, help="Path to the config file.")
    parser.add_argument("--ext", dest="data_file_ext", default='ecsv',
                        type=str, help="Extension of data files.")

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

    # Load data here
    filename = path.join(args.data_path, '{0}.{1}'.format(args.apogee_id,
                                                          args.data_file_ext))
    data_tbl = QTable.read(filename)
    data = RVData(t=data_tbl['time'], rv=data_tbl['rv'],
                  stddev=data_tbl['rv_err'])

    # parse config file
    with open(args.config_file, 'r') as f:
        config = yaml.load(f.read())
        config['config_file'] = args.config_file
    joker_params = config_to_jokerparams(config)

    pool_kwargs = dict(mpi=args.mpi, processes=args.n_procs)
    pool = choose_pool(**pool_kwargs)

    main(data_path=args.data_path, pool=pool,
         apogee_id=args.apogee_id,
         config=config, overwrite=args.overwrite,
         data_file_ext=args.data_file_ext)

    pool.close()
    sys.exit(0)
