# Standard library
from os.path import abspath, expanduser

# Third-party
from astropy.io import fits
from astropy.table import Table
import numpy as np

# Project
from twoface.util import Timer
from twoface.log import log as logger
from twoface.db.query_helpers import paged_query

from apogeebh.db.connect import db_connect, Base
from apogeebh.db.model import AllStar, AllVisit, Status

__all__ = ['initialize_db']


def tblrow_to_dbrow(tblrow, colnames, varchar_cols=[]):
    row_data = dict()
    for k in colnames:
        if k in varchar_cols:
            row_data[k.lower()] = tblrow[k].strip()

        # special case the array columns
        elif k.lower().startswith('fparam') and 'cov' not in k.lower():
            i = int(k[-1])
            row_data[k.lower()] = tblrow[k[:-1]][i]

        elif k.lower().startswith('fparam') and 'cov' in k.lower():
            i = int(k[-2])
            j = int(k[-1])
            row_data[k.lower()] = tblrow[k[:-2]][i,j]

        # foreign keys:
        elif k.lower() in ['allstar_id']:
            continue

        else:
            row_data[k.lower()] = tblrow[k]

    return row_data


def initialize_db(allVisit_file, allStar_file, database_path,
                  drop_all=False, batch_size=4096):
    """Initialize the database given FITS filenames for the APOGEE data.

    Parameters
    ----------
    allVisit_file : str
        Full path to APOGEE allVisit file.
    allStar_file : str
        Full path to APOGEE allStar file.
    database_file : str
        Filename (not path) of database file in cache path.
    drop_all : bool (optional)
        Drop all existing tables and re-create the database.
    batch_size : int (optional)
        How many rows to create before committing.
    """

    norm = lambda x: abspath(expanduser(x))
    allvisit = fits.getdata(norm(allVisit_file))
    allstar = fits.getdata(norm(allStar_file))

    # STAR_BAD
    ASPCAP_skip_bitmask = np.sum(2 ** np.array([23]))

    # VERY_BRIGHT_NEIGHBOR, LOW_SNR, SUSPECT_RV_COMBINATION, SUSPECT_BROAD_LINES
    STAR_skip_bitmask = np.sum(2 ** np.array([3, 4, 16, 17]))

    # First filter allStar flags
    mask = ((allstar['NVISITS'] >= 4) &
            ((allstar['ASPCAPFLAG'] & ASPCAP_skip_bitmask) == 0) &
            ((allstar['STARFLAG'] & STAR_skip_bitmask) == 0)
           )
    stars = allstar[mask]
    visits = allvisit[np.isin(allvisit['APOGEE_ID'], stars['APOGEE_ID'])]

    # Next filter allVisit flags
    mask = (((visits['STARFLAG'] & STAR_skip_bitmask) == 0) &
            (visits['VRELERR'] < 100.) &
            np.isfinite(visits['VHELIO']) & np.isfinite(visits['VRELERR']) &
            (visits['VHELIO'] > -999)
           )
    visits = visits[mask]
    v_apogee_ids, counts = np.unique(visits['APOGEE_ID'], return_counts=True)
    stars = stars[np.isin(stars['APOGEE_ID'], v_apogee_ids[counts >= 4])]
    visits = visits[np.isin(visits['APOGEE_ID'], stars['APOGEE_ID'])]

    # uniquify the stars
    _, idx = np.unique(stars['APOGEE_ID'], return_index=True)
    allstar_tbl = Table(stars[idx])
    allvisit_tbl = Table(visits)

    Session, engine = db_connect(database_path, ensure_db_exists=True)
    logger.debug("Connected to database at '{}'".format(database_path))

    if drop_all:
        # this is the magic that creates the tables based on the definitions in
        # twoface/db/model.py
        Base.metadata.drop_all()
        Base.metadata.create_all()

    session = Session()

    logger.debug("Loading allStar, allVisit tables...")

    # Figure out what data we need to pull out of the FITS files based on what
    # columns exist in the (empty) database
    allstar_skip = ['ID']
    allstar_colnames = []
    allstar_varchar = []
    for x in AllStar.__table__.columns:
        col = str(x).split('.')[1].upper()
        if col in allstar_skip:
            continue

        if str(x.type) == 'VARCHAR':
            allstar_varchar.append(col)

        allstar_colnames.append(col)

    allvisit_skip = ['ID']
    allvisit_colnames = []
    allvisit_varchar = []
    for x in AllVisit.__table__.columns:
        col = str(x).split('.')[1].upper()
        if col in allvisit_skip:
            continue

        if str(x.type) == 'VARCHAR':
            allvisit_varchar.append(col)

        allvisit_colnames.append(col)

    # --------------------------------------------------------------------------
    # First load the status table:
    #
    if session.query(Status).count() == 0:
        logger.debug("Populating Status table...")
        statuses = list()
        statuses.append(Status(id=0, message='untouched'))
        statuses.append(Status(id=1, message='needs more prior samples'))
        statuses.append(Status(id=2, message='needs mcmc'))
        statuses.append(Status(id=3, message='error'))
        statuses.append(Status(id=4, message='completed'))

        session.add_all(statuses)
        session.commit()
        logger.debug("...done")

    # --------------------------------------------------------------------------
    # Load the AllStar table:
    #
    logger.info("Loading AllStar table")

    # What APOGEE_ID's are already loaded?
    all_ap_ids = np.array([x.strip() for x in allstar_tbl['APOGEE_ID']])
    loaded_ap_ids = [x[0] for x in session.query(AllStar.apogee_id).all()]
    mask = np.logical_not(np.isin(all_ap_ids, loaded_ap_ids))
    logger.debug("{0} stars already loaded".format(len(loaded_ap_ids)))
    logger.debug("{0} stars left to load".format(mask.sum()))

    stars = []
    with Timer() as t:
        i = 0
        for row in allstar_tbl[mask]: # Load every star
            row_data = tblrow_to_dbrow(row, allstar_colnames, allstar_varchar)

            # create a new object for this row
            star = AllStar(**row_data)
            stars.append(star)
            logger.log(1, 'Adding star {0} to database'.format(star))

            if i % batch_size == 0 and i > 0:
                session.add_all(stars)
                session.commit()
                logger.debug("Loaded batch {0} ({1:.2f} seconds)"
                             .format(i, t.elapsed()))
                t.reset()
                stars = []

            i += 1

    if len(stars) > 0:
        session.add_all(stars)
        session.commit()

    # --------------------------------------------------------------------------
    # Load the AllVisit table:
    #
    logger.info("Loading AllVisit table")

    # What VISIT_ID's are already loaded?
    all_vis_ids = np.array([x.strip() for x in allvisit_tbl['VISIT_ID']])
    loaded_vis_ids = [x[0] for x in session.query(AllVisit.visit_id).all()]
    mask = np.logical_not(np.isin(all_vis_ids, loaded_vis_ids))
    logger.debug("{0} visits already loaded".format(len(loaded_vis_ids)))
    logger.debug("{0} visits left to load".format(mask.sum()))

    visits = []
    with Timer() as t:
        i = 0
        for row in allvisit_tbl[mask]: # Load every visit
            row_data = tblrow_to_dbrow(row, allvisit_colnames, allvisit_varchar)

            # create a new object for this row
            visit = AllVisit(**row_data)
            visits.append(visit)
            logger.log(1, 'Adding visit {0} to database'.format(visit))

            if i % batch_size == 0 and i > 0:
                session.add_all(visits)
                session.commit()
                logger.debug("Loaded batch {0} ({1:.2f} seconds)"
                             .format(i, t.elapsed()))
                t.reset()
                visits = []

            i += 1

    if len(visits) > 0:
        session.add_all(visits)
        session.commit()

    # --------------------------------------------------------------------------
    # Now associate rows in AllStar with rows in AllVisit
    logger.info("Linking AllVisit and AllStar tables")

    q = session.query(AllStar).order_by(AllStar.id)

    for i, sub_q in enumerate(paged_query(q, page_size=batch_size)):
        for star in sub_q:
            if len(star.visits) > 0:
                continue

            visits = session.query(AllVisit).filter(
                AllVisit.apogee_id == star.apogee_id).all()

            if len(visits) == 0:
                logger.warn("Visits not found for star {0}".format(star))
                continue

            logger.log(1, 'Attaching {0} visits to star {1}'
                       .format(len(visits), star))

            star.visits = visits

        logger.debug("Committing batch {0}".format(i))
        session.commit()

    session.commit()
    session.close()
