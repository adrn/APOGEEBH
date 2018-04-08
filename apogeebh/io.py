# Third-party
from astropy.table import QTable


def load_data_ecsv(fn):
    tbl = QTable.read(fn)
    data = RVData(t=tbl['time'],
                  rv=tbl['rv'],
                  stddev=tbl['rv_err'])
    return data, tbl['source']
