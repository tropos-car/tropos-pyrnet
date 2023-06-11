# AUTOGENERATED! DO NOT EDIT! File to edit: ../../nbs/pyrnet/data.ipynb.

# %% auto 0
__all__ = ['pyrnet_version', 'logger', 'update_coverage_meta', 'stretch_resolution', 'merge_ds', 'to_netcdf', 'get_config',
           'get_cfmeta', 'to_l1a', 'to_l1b']

# %% ../../nbs/pyrnet/data.ipynb 2
import os
import numpy as np
import pandas as pd
import xarray as xr
import logging
from toolz import assoc_in
import pkg_resources as pkg_res
import warnings

from trosat import sunpos as sp

import pyrnet as pyrnet_main
pyrnet_version = pyrnet_main.__version__
from . import pyrnet
from . import utils as pyrutils
from . import logger as pyrlogger
from . import reports as pyrreports

# logging setup
logging.basicConfig(
    filename='pyrnet.log',
    encoding='utf-8',
    level=logging.DEBUG,
    format='%(asctime)s %(name)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)

# %% ../../nbs/pyrnet/data.ipynb 5
def update_coverage_meta(ds,timevar='time'):
    """Update global attributes related to geospatial and time coverage
    """
    duration = ds[timevar].values[-1] - ds[timevar].values[0]
    resolution = np.mean(np.diff(ds[timevar].values))
    now = pd.to_datetime(np.datetime64("now"))
    gattrs = {
        'date_created': now.isoformat(),
        'geospatial_lat_min': np.nanmin(ds.lat.values),
        'geospatial_lat_max': np.nanmax(ds.lat.values),
        'geospatial_lat_units': 'degN',
        'geospatial_lon_min': np.nanmin(ds.lon.values),
        'geospatial_lon_max': np.nanmax(ds.lon.values),
        'geospatial_lon_units': 'degE',
        'time_coverage_start': pd.to_datetime(ds[timevar].values[0]).isoformat(),
        'time_coverage_end': pd.to_datetime(ds[timevar].values[-1]).isoformat(),
        'time_coverage_duration': pd.to_timedelta(duration).isoformat(),
        'time_coverage_resolution': pd.to_timedelta(resolution).isoformat(),
    }
    ds.attrs.update(gattrs)
    return ds


# %% ../../nbs/pyrnet/data.ipynb 6
def stretch_resolution(ds: xr.Dataset) -> xr.Dataset:
    """ Stretch variable resolution to full integer size,
    to not lose resolution after averaging ADC count data."""
    for var in ds:
        dtype = ds[var].encoding['dtype']
        valid_range = ds[var].valid_range
        int_limit = np.iinfo(dtype).max
        scale_factor = ds[var].encoding['scale_factor']
        scale_factor_mod = int((int_limit-1)/valid_range[1])
        ds[var].encoding.update({
            "scale_factor": scale_factor * scale_factor_mod,
            "_FillValue": int_limit,
            "valid_range": valid_range * scale_factor_mod
        })
    return ds

# %% ../../nbs/pyrnet/data.ipynb 8
def merge_ds(ds1,ds2):
    """Merge two datasets along the time dimension.
    """
    if ds1.time.equals(ds2.time):
        logging.info("Overwrite existing file.")
        return ds2
    logging.info("Merge with existing file.")

    ## overwrite non time dependent variables
    overwrite_vars = [ v for v in ds1 if "time" not in ds1[v].dims]

    ## merge both datasets
    ds_new=ds1.merge(ds2,compat='no_conflicts',overwrite_vars=overwrite_vars)

    # add global coverage attributes
    ds_new.attrs.update({'merged':1})
    return ds_new

# %% ../../nbs/pyrnet/data.ipynb 9
def to_netcdf(ds,fname):
    """xarray to netcdf, but merge if exist
    """
    # merge if necessary
    if os.path.exists(fname):
        ds1 = xr.open_dataset(fname)
        ds = merge_ds(ds1,ds)
        os.remove(fname)

    # save to netCDF4
    ds = update_coverage_meta(ds, timevar="time")
    ds.to_netcdf(fname,
                 encoding={'time':{'dtype':'float64'}}) # for OpenDAP 2 compatibility

# %% ../../nbs/pyrnet/data.ipynb 11
def get_config(config: dict|None = None) -> dict:
    """Read default config and merge with input config
    """

    fn_config = pkg_res.resource_filename("pyrnet", "share/pyrnet_config.json")
    default_config = pyrutils.read_json(fn_config)
    if config is None:
        config = default_config
    config = {**default_config, **config}

    # add default files
    cfiles = {
        "file_cfmeta": "share/pyrnet_cfmeta.json",
        "file_calibration": "share/pyrnet_calibration.json",
        "file_mapping": "share/pyrnet_station_map.json",
        "file_gti_angles": "share/pyrnet_gti_angles.json",
        "file_site": "share/pyrnet_sites.json",
    }
    for fn in cfiles:
        if config[fn] is None:
            config[fn] =  pkg_res.resource_filename("pyrnet", cfiles[fn])
    return config

def get_cfmeta(config: dict|None = None) -> dict:
    """Read global and variable attributes and encoding from cfmeta.json
    """
    config= get_config(config)
    # parse the json file
    cfdict = pyrutils.read_json(config["file_cfmeta"])
    # get global attributes:
    gattrs = cfdict['attributes']
    # apply config
    gattrs = {k:v.format_map(config) for k,v in gattrs.items()}
    # get variable attributes
    d = pyrutils.get_var_attrs(cfdict)
    # split encoding attributes
    vattrs, vencode = pyrutils.get_attrs_enc(d)
    return gattrs ,vattrs, vencode


# %% ../../nbs/pyrnet/data.ipynb 17
def to_l1a(
        fname : str,
        *,
        station: int,
        report: dict,
        date_of_measure : np.datetime64 = np.datetime64("now"),
        config: dict|None = None,
        global_attrs: dict|None = None
) -> xr.Dataset|None:
    """
    Read logger raw file and parse it to xarray Dataset. Thereby, attributes and names are defined via cfmeta.json file and sun position values are calculated and added.

    Parameters
    ----------
    fname: str
        Path and filename of the raw logger file.
    station: int
        PyrNet station box number.
    report: dict
        Parsed maintenance report, see reports.ipynb
    bins: int
        Number of desired bins per day. The default is 86400, which result in
        mean values of 1 second steps per day. Maximum resolution is 86400000.
    date_of_measure: float, datetime or datetime64
        A rough date of measure  to account for GPS week rollover. If measured in 2019, day resolution is recommended, before 2019 annual resolution, 2020 onwards not required. If float, interpreted as Julian day from 2000-01-01T12:00. the default is np.datetime64("now").
    config: dict
        Stores processing specific configuration.
            * cfjson -> path to cfmeta.json, the default is "../share/pyrnet_cfmeta.json"
            * stripminutes -> number of minutes to be stripped from the data at start and end,
                the default is 5.
    global_attrs: dict
        Additional global attributes for the Dataset. (Overrides cfmeta.json attributes)
    Returns
    -------
    xarray.Dataset
        Raw Logger data for one measurement periode.
    """
    config = get_config(config)
    gattrs, vattrs, vencode = get_cfmeta(config)

    if global_attrs is not None:
        gattrs.update(global_attrs)

    date_of_measure = pyrutils.to_datetime64(date_of_measure)

    # 1. Parse raw file
    rec_adc, rec_gprmc = pyrlogger.read_records(fname=fname, date_of_measure=date_of_measure)

    if type(rec_adc)==bool or len(rec_gprmc.time)<3:
        logger.debug("Failed to load the data from the file, because of not enough stable GPS data, or file is empty.")
        return None

    # Get ADC time
    adctime = pyrlogger.get_adc_time(rec_adc)

    # ADC to Volts
    # Drop time and internal battery sensor output (columns 0 and 1)
    adc_volts = 3.3 * rec_adc[:,2:] / 1023.

    # 2. Get Logbook maintenance quality flags
    key = f"{station:03d}"
    if isinstance(report, pd.DataFrame):
        logger.info(f"Parsing report at date {rec_gprmc.time[-1]}")
        report = pyrreports.parse_report(report,
                                date_of_maintenance=rec_gprmc.time[-1])

    if key not in report:
        logger.warning(f"No report for station {station} available.")
        warnings.warn(f"No report for station {station} available.")
        qc_main = pyrreports.get_qcflag(4,3)
        qc_extra = pyrreports.get_qcflag(4,3)
        vattrs = assoc_in(vattrs, ["ghi_qc","note_general"], "No maintenance report!")
        vattrs = assoc_in(vattrs, ["gti_qc","note_general"], "No maintenance report!")
    else:
        qc_main = pyrreports.get_qcflag(
            qc_clean=report[key]['clean'],
            qc_level=report[key]['align']
        )
        qc_extra = pyrreports.get_qcflag(
            qc_clean=report[key]['clean2'],
            qc_level=report[key]['align2']
        )
        # add qc notes
        vattrs = assoc_in(vattrs, ["ghi_qc","note_general"], report[key]["note_general"])
        vattrs = assoc_in(vattrs, ["gti_qc","note_general"], report[key]["note_general"])
        vattrs = assoc_in(vattrs, ["ghi_qc","note_clean"], report[key]["note_clean"])
        vattrs = assoc_in(vattrs, ["gti_qc","note_clean"], report[key]["note_clean2"])
        vattrs = assoc_in(vattrs, ["ghi_qc","note_level"], report[key]["note_align"])
        vattrs = assoc_in(vattrs, ["gti_qc","note_level"], report[key]["note_align2"])

    # 3. Add global meta data
    now = pd.to_datetime(np.datetime64("now"))
    gattrs.update({
        'processing_level': 'l1a',
        'product_version': pyrnet_version,
        'history': f'{now.isoformat()}: Generated level l1a  by pyrnet version {pyrnet_version}; ',
    })
    # add site information
    if config['sites'] is not None:
        sites = pyrutils.read_json(config['file_site'])[config['sites']]
        if key in sites:
            gattrs.update({ "site" : sites[key]})



    # add gti angles
    # default horizontal
    vattrs = assoc_in(vattrs, ["gti","hangle"], 0)
    vattrs = assoc_in(vattrs, ["gti","vangle"], 0)
    # update with angles from mapping file
    if config['gti_angles'] is not None:
        gti_angles = pyrutils.read_json(config['file_gti_angles'])[config['gti_angles']]
        if key in gti_angles:
            vattrs = assoc_in(vattrs, ["gti","hangle"], gti_angles[key][0])
            vattrs = assoc_in(vattrs, ["gti","vangle"], gti_angles[key][1])

    if adc_volts.shape[1]<5: # gti data is not available
        adc_volts = np.concatenate((adc_volts,-1*np.ones(adc_volts.shape[0])[:,None]),axis=1)

    # 8. Make xarray Dataset
    ds = xr.Dataset(
        data_vars={
            "ghi": (("adctime","station"), adc_volts[:,2][:,None] / 300.), # [V]
            "gti": (("adctime","station"), adc_volts[:,4][:,None] / 300.), # [V]
            "ta": (("adctime","station"), 253.15 + 20.*2.*adc_volts[:,0][:,None]), # [K]
            "rh": (("adctime","station"), 0.2*2.*adc_volts[:,1][:,None]), # [-]
            "battery_voltage": (("adctime","station"), 2.*adc_volts[:,3][:,None]), # [V]
            "lat": (("gpstime","station"), rec_gprmc.lat[:,None]), # [degN]
            "lon": (("gpstime","station"), rec_gprmc.lon[:,None]), # [degE]
            "ghi_qc": ("station", [qc_main]),
            "gti_qc": ("station", [qc_extra]),
            "iadc": (("gpstime", "station"), rec_gprmc.iadc[:,None])
        },
        coords={
            "adctime": ("adctime", adctime.astype('timedelta64[ns]')),
            "gpstime": ("gpstime", rec_gprmc.time.astype('datetime64[ns]')),
            "station": ("station", [station]),
        },
        attrs=gattrs
    )

    # drop ocurance of douplicate gps values
    ds = ds.drop_duplicates("gpstime")

    # add global coverage attributes
    ds = update_coverage_meta(ds, timevar="gpstime")

    # add attributes to Dataset
    for k,v in vattrs.items():
        if k not in ds.keys():
            continue
        ds[k].attrs = v

    # add encoding to Dataset
    for k,v in vencode.items():
        if k not in ds.keys():
            continue
        ds[k].encoding = v

    ds["gpstime"].encoding.update({
        "dtype": 'f8',
        "units": f"seconds since {np.datetime_as_string(ds.gpstime.data[0], unit='D')}T00:00Z",
    })
    ds["adctime"].encoding.update({
        "units": f"milliseconds",
    })

    return ds

# %% ../../nbs/pyrnet/data.ipynb 44
def to_l1b(
        fname: str,
        *,
        config: dict | None = None,
        global_attrs: dict | None = None
) -> xr.Dataset:

    config = get_config(config)
    gattrs, vattrs, vencode = get_cfmeta(config)

    if global_attrs is not None:
        gattrs.update(global_attrs)

    # 1. Load l1a data
    ds_l1a = xr.open_dataset(fname)
    # check correct file
    if ds_l1a.processing_level != "l1a":
        raise ValueError(f"{fname} is not a l1a file.")

    # 2. Sync GPS to ADC time
    adctime = pyrlogger.sync_adc_time(
        adctime = ds_l1a.adctime.values,
        gpstime = ds_l1a.gpstime.values,
        iadc = ds_l1a.iadc.squeeze().values.astype(int)
    )

    # 3. Create new dataset (l1b)
    ds_l1b = ds_l1a.drop_dims('gpstime')
    ds_l1b = ds_l1b.drop_vars(['ghi_qc','gti_qc']) # keep only time dependend variables
    ds_l1b = ds_l1b.assign({'time': ('adctime', adctime)})
    ds_l1b = ds_l1b.swap_dims({"adctime":"time"})
    ds_l1b = ds_l1b.drop_vars("adctime")

    ds_l1b["time"].encoding.update({
        "dtype": 'float64',
        "units": f"seconds since {np.datetime_as_string(ds_l1b.time.data[0], unit='D')}T00:00Z",
    })
    logger.info(f"Dataset time coverage before strip: {ds_l1b.time.values[0]} - {ds_l1b.time.values[-1]}")

    # 4. Drop first and last <stripminutes> minutes of data to avoid bad data due to maintenance
    stripminutes = np.timedelta64(int(config['stripminutes']), 'm')
    ds_l1b = ds_l1b.isel(time=ds_l1b.time>ds_l1b.time.values[0] + stripminutes)
    ds_l1b = ds_l1b.isel(time=ds_l1b.time<ds_l1b.time.values[-1] - stripminutes)
    logger.info(f"Dataset time coverage after strip: {ds_l1b.time.values[0]} - {ds_l1b.time.values[-1]}")


    # 5. resample to desired resolution
    ds_l1b = pyrlogger.resample_mean(ds_l1b, freq=config['l1bfreq'])
    # stretch valid range to not lose resolution due to averaging
    ds_l1b = stretch_resolution(ds_l1b)

    # 6. Interpolate GPS coordinates to l1b time
    ds_gps = ds_l1a.drop_dims("adctime")
    ds_gps = ds_gps.drop_vars(['iadc'])

    # Decide whether geo coordinates should be averaged or not
    if config['average_latlon']:
        ds_gps = ds_gps.mean('gpstime', skipna=True, keep_attrs=True)
    else:
        ds_gps = ds_gps.interp(gpstime=ds_l1b.time)
        ds_gps = ds_gps.drop_vars("gpstime")

    ds_l1b = xr.merge((ds_l1b,ds_gps))

    # 7. Calc and add sun position
    szen, sazi = sp.sun_angles(
        time=ds_l1b.time.values[:,None], # line up with coordinates to keep dependence on time only
        lat=ds_l1b.lat.values,
        lon=ds_l1b.lon.values
    )
    szen  = szen.squeeze()
    sazi = sazi.squeeze()

    esd = np.mean(sp.earth_sun_distance(ds_l1b.time.values))

    ds_l1b = ds_l1b.assign(
        {
            "szen": (("time", "station"), szen[:,None]),
            "sazi": (("time", "station"), sazi[:,None]),
            "esd": ("station", [esd])
        }
    )
    # update attributes and encoding
    for key in ['szen', 'sazi','esd']:
        ds_l1b[key].attrs.update(vattrs[key])
        ds_l1b[key].encoding.update(vencode[key])


    # 8. rad flux calibration
    box = ds_l1b.station.values[0]
    boxnumber, serial, cfac = pyrnet.meta_lookup(
        ds_l1b.time.values[0],
        box=box,
        cfile=config['file_calibration'],
        mapfile=config['file_mapping'],
    )
    logger.info(f"Meta Lookup:")
    logger.info(f">> Box={box}")
    logger.info(f">> serial(s)={serial}")
    logger.info(f">> calibration factor(s)={cfac}")


    # calibrate radiation flux with gain=300
    for i, radflx in enumerate(config['radflux_varname']):
        if cfac[i] is None:
            # drop if calibration/instrument don't exist (probably secondary pyranometer).
            ds_l1b = ds_l1b.drop_vars([var for var in ds_l1b if radflx in var])
            continue
        ds_l1b[radflx].values = ds_l1b[radflx].values*1e6/(cfac[i]) # V -> W m-2
        ds_l1b[radflx].attrs['units'] = "W m-2",
        ds_l1b[radflx].attrs.update({
            "units": "W m-2",
            "serial": serial[i],
            "calibration_factor": cfac[i]
        })
        ds_l1b[radflx].encoding.update({
            'scale_factor': ds_l1b[radflx].encoding['scale_factor']*1e6/(cfac[i])
        })


    # add global coverage attributes
    ds_l1b = update_coverage_meta(ds_l1b, timevar="time")

    ds_l1b.attrs["processing_level"] = 'l1b'
    now = pd.to_datetime(np.datetime64("now"))
    ds_l1b.attrs["history"] = ds_l1b.history + f"{now.isoformat()}: Generated level l1b  by pyrnet version {pyrnet_version}; "

    return ds_l1b
