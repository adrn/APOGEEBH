# Project
from twoface.log import log as logger
from apogeebh.db.init import initialize_db

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="Initialize the TwoFace project "
                                        "database.")

    vq_group = parser.add_mutually_exclusive_group()
    vq_group.add_argument('-v', '--verbose', action='count', default=0,
                          dest='verbosity')
    vq_group.add_argument('-q', '--quiet', action='count', default=0,
                          dest='quietness')

    parser.add_argument('-o', '--overwrite', action='store_true', default=False,
                        dest='overwrite')

    parser.add_argument("--allstar", dest="allStar_file", required=True,
                        type=str, help="Path to APOGEE allStar FITS file.")
    parser.add_argument("--allvisit", dest="allVisit_file", required=True,
                        type=str, help="Path to APOGEE allVisit FITS file.")
    parser.add_argument("--db", dest="db_file", required=True,
                        type=str, help="Path to database file.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbosity != 0:
        if args.verbosity == 1:
            logger.setLevel(logging.DEBUG)
        else: # anything >= 2
            logger.setLevel(1)

    elif args.quietness != 0:
        if args.quietness == 1:
            logger.setLevel(logging.WARNING)
        else: # anything >= 2
            logger.setLevel(logging.ERROR)

    else: # default
        logger.setLevel(logging.INFO)

    initialize_db(args.allVisit_file, args.allStar_file, args.db_file,
                  drop_all=args.overwrite)
