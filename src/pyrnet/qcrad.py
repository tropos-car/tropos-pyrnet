# AUTOGENERATED! DO NOT EDIT! File to edit: ../../nbs/pyrnet/qcrad.ipynb.

# %% auto 0
__all__ = ['logger', 'CONSTANTS', 'QCCode', 'init_qc_flag']

# %% ../../nbs/pyrnet/qcrad.ipynb 2
import xarray as xr
import numpy as np
import trosat.sunpos as sp
import logging

import pyrnet.data


# logging setup
logging.basicConfig(
    filename='pyrnet.log',
    encoding='utf-8',
    level=logging.DEBUG,
    format='%(asctime)s %(name)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)

# %% ../../nbs/pyrnet/qcrad.ipynb 3
class CONSTANTS:
    S0 = 1367  # W m-2
    k = 5.67*1e-8

# %% ../../nbs/pyrnet/qcrad.ipynb 5
class QCCode:
    """ BSRN quality codes
    https://wiki.pangaea.de/wiki/BSRN_Toolbox#Quality_Check
    """
    below_physical = 2**0
    above_phyiscal = 2**1
    below_rare = 2**2
    above_rare = 2**3
    compare_to_low = 2**4
    compare_to_high = 2**5

# %% ../../nbs/pyrnet/qcrad.ipynb 9
def init_qc_flag(ds, var):
    qc_bits = [2**i for i in range(6)]
    ds[f"qc_flag_{var}"] = ds[var].copy()
    ds[f"qc_flag_{var}"].values *= 0
    ds[f"qc_flag_{var}"].attrs.update({
        "standard_name": "quality_flag",
        "units": "-",
        "ancillary_variables": var,
        "valid_range": [0, np.sum(qc_bits)],
        "flag_masks": qc_bits,
        "flag_values": qc_bits,
        "flag_meanings": str(
            "below_physical_limit" + " " +
            "above_physical_limit" + " " +
            "below_rare_limit" + " " +
            "above_rare_limit" + " " +
            "comparison_to_low" + " " +
            "comparison_to_high"
        )
    })
    ds[f"qc_flag_{var}"].encoding.update({
        "dtype": "u1",
        "_FillValue": 0
    })
    return ds
