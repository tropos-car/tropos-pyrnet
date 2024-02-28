# AUTOGENERATED! DO NOT EDIT! File to edit: ../../nbs/pyrnet/qcrad.ipynb.

# %% auto 0
__all__ = ['logger', 'CONSTANTS', 'QCCode', 'init_qc_flag', 'add_qc_flags']

# %% ../../nbs/pyrnet/qcrad.ipynb 2
import xarray as xr
import numpy as np
import warnings
import logging

import pyrnet.data
import pyrnet.utils


# logging setup
logging.basicConfig(
    filename='pyrnet.log',
    encoding='utf-8',
    level=logging.DEBUG,
    format='%(asctime)s %(name)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)

# %% ../../nbs/pyrnet/qcrad.ipynb 4
class CONSTANTS:
    S0 = 1367  # W m-2
    k = 5.67*1e-8

# %% ../../nbs/pyrnet/qcrad.ipynb 6
class QCCode:
    """ BSRN quality codes
    https://wiki.pangaea.de/wiki/BSRN_Toolbox#Quality_Check
    """
    below_physical = 2**0 # 1
    above_phyiscal = 2**1 # 2
    below_rare = 2**2 # 4
    above_rare = 2**3 # 8
    compare_to_low = 2**4 # 16
    compare_to_high = 2**5 # 32
    quality_control_failed = 2**6 # 64

# %% ../../nbs/pyrnet/qcrad.ipynb 10
def init_qc_flag(ds, var):
    qc_bits = [2**i for i in range(7)]
    ds[f"qc_flag_{var}"] = ds[var].copy()
    ds[f"qc_flag_{var}"].values = np.zeros(ds[var].shape).astype(np.ubyte)
    attrs = ds[f"qc_flag_{var}"].attrs
    attrs.update({
        "standard_name": "quality_flag",
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
            "comparison_to_high" + " "+
            "quality_control_failed"
        )
    })
    attrs.pop("units", None)
    attrs.pop("long_name", None)
    ds[f"qc_flag_{var}"].attrs.update(attrs)
    ds[f"qc_flag_{var}"].encoding.update({
        "dtype": "u1",
        "_FillValue": 255,
        "zlib": True
    })
    return ds

# %% ../../nbs/pyrnet/qcrad.ipynb 18
def add_qc_flags(ds, vars):
    """
    Add quality flags to flux variables in the dataset.
    
    Parameters
    ----------
    ds: xr.Dataset
        Dataset with flux variables, dimensions ('time','station'). Also, solar zenith (szen) and azimuth (sazi) angles are required.
        Works for pyrnet l1b and l1b_network files.
    vars: list
        List of flux variable names in ds.

    Returns
    -------
    xr.Dataset
        The input dataset, but with additional 'qc_flag_<fluxvar>' variables.
    """
    # keep only available variables
    vars = [ var for var in vars if var in ds ]
    
    # init qc flags
    for var in vars:
        ds = init_qc_flag(ds, var) 
    
    # ancillary variables
    szen = ds.szen.values
    mu0 = np.cos(np.deg2rad(szen))
    mu0[mu0 < 0] = 0 #  exclude night
    esd = ds.esd.values
    Sa = CONSTANTS.S0 / esd**2
    
    # prepare subsample dataset, to be resampled to 30min mean for comparison checks
    dsr = ds.copy()
    dsr = dsr.drop_vars([var for var in dsr if var not in vars])
    
    # do physical and extreme limit tests
    for var in vars:
        is_tilted = pyrnet.utils.check_tilted(ds[var])
        values = ds[var].values.copy()
        
        # apply correction if possible
        if np.any(is_tilted):
            vangle = pyrnet.utils.make_iter(ds[var].attrs["vangle"])
            hangle = pyrnet.utils.make_iter(ds[var].attrs["hangle"])
            cfac = pyrnet.utils.tilt_correction_factor(
                dp = vangle,
                dy = hangle,
                szen=ds.szen.values,
                sazi=ds.sazi.values
            )
            mask = is_tilted[None,:] * np.isnan(cfac)
            ds[f"qc_flag_{var}"].values[mask] += QCCode.quality_control_failed
            apply_correction = is_tilted[None,:] * ~np.isnan(cfac)
            values[apply_correction] *= cfac[apply_correction]
        
        # update subsample dataset
        dsr[var].values = values
        
        # physical minimum
        mask = values < -4
        ds[f"qc_flag_{var}"].values[mask] += QCCode.below_physical
        # physical maximum
        mask = values > ((Sa * 1.5 * mu0 ** 1.2) + 100)
        ds[f"qc_flag_{var}"].values[mask] += QCCode.above_phyiscal
        # rare limit minimum
        mask = values < -2
        ds[f"qc_flag_{var}"].values[mask] += QCCode.below_rare
        # rare limit maximum
        mask = values > ((Sa * 1.2 * mu0 ** 1.2) + 50)
        ds[f"qc_flag_{var}"].values[mask] += QCCode.above_rare
        
        
    # compare all sensors from network, or single station
    dsr = dsr.resample(time="30min",skipna=True).mean()
    thres_low = np.ones(ds.time.size)*0.9
    thres_high = np.ones(ds.time.size)*1.1
    thres_low[ds.szen.mean("station")>75] = 0.85
    thres_low[ds.szen.mean("station")>75] = 1.15
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        all_values_mean = np.nanmean(np.concatenate([dsr[var].values for var in dsr],axis=1),axis=1)
        
    for var in vars:
        is_tilted = pyrnet.utils.check_tilted(ds[var])
        ratio = np.ones(dsr[var].shape)
        if np.any(is_tilted):
            meanvalues =  all_values_mean
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                meanvalues = np.nanmean(dsr[var].values, axis=1)
        ratio[meanvalues>50] = dsr[var].values[meanvalues>50] / meanvalues[meanvalues>50][:,None]
        
        # reindex ratio to original resolution
        dsr = dsr.assign({"ratio": (("time","station"), ratio)})
        ratio = dsr.ratio.reindex_like(ds, method='nearest').values
        dsr = dsr.drop_vars(["ratio"])
        # comparison to low
        mask = ratio < thres_low[:,None]
        ds[f"qc_flag_{var}"].values[mask] += QCCode.compare_to_low
        
        # comparison to high
        mask = ratio > thres_high[:,None]
        ds[f"qc_flag_{var}"].values[mask] += QCCode.compare_to_high
    
    return ds
