# AUTOGENERATED! DO NOT EDIT! File to edit: ../../nbs/pyrnet/data.ipynb.

# %% auto 0
__all__ = ['pyrnet_version', 'logger', 'update_coverage_meta', 'stretch_resolution', 'merge_ds', 'to_netcdf', 'resample',
           'get_config', 'get_sensor_config', 'get_cfmeta', 'calc_encoding', 'add_encoding', 'to_l1a', 'to_l1b']

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
        if "scale_factor" not in ds[var].encoding:
            continue
        if "valid_range" not in ds[var].attrs:
            continue
        dtype = ds[var].encoding['dtype']
        valid_range = ds[var].valid_range
        int_limit = np.iinfo(dtype).max
        scale_factor = ds[var].encoding['scale_factor']
        scale_factor_mod = int((int_limit-1)/valid_range[1])
        ds[var].encoding.update({
            "scale_factor": scale_factor / scale_factor_mod,
            "_FillValue": int_limit,
        })
        ds[var].attrs.update({
            "valid_range": valid_range * scale_factor_mod
        })
    return ds

# %% ../../nbs/pyrnet/data.ipynb 7
def merge_ds(ds1,ds2,timevar="time"):
    """Merge two datasets along the time dimension.
    """
    if ds1[timevar].equals(ds2[timevar]):
        logging.info("Overwrite existing file.")
        return ds2
    logging.info("Merge with existing file.")

    ## overwrite non time dependent variables
    overwrite_vars = [ v for v in ds1 if timevar not in ds1[v].dims ]

    ## merge both datasets
    ds_new=ds1.merge(ds2,
                     compat='no_conflicts',
                     overwrite_vars=overwrite_vars)

    # add global coverage attributes
    ds_new.attrs.update({'merged':1})

    # add encoding again
    ds_new = add_encoding(ds_new)
    return ds_new

# %% ../../nbs/pyrnet/data.ipynb 8
def to_netcdf(ds,fname, timevar="time"):
    """xarray to netcdf, but merge if exist
    """
    # merge if necessary
    if os.path.exists(fname):
        ds1 = xr.open_dataset(fname)
        ds = merge_ds(ds1,ds, timevar=timevar)
        ds1.close()
        os.remove(fname)

    # save to netCDF4
    ds = update_coverage_meta(ds, timevar=timevar)
    ds.to_netcdf(fname,
                 encoding={timevar:{'dtype':'float64'}}) # for OpenDAP 2 compatibility

# %% ../../nbs/pyrnet/data.ipynb 9
def resample(ds, freq, methods='mean', kwargs={}):
    """ Resample xarray dataset using pandas for speed.
    https://github.com/pydata/xarray/issues/4498#issuecomment-706688398
    """
    if isinstance(methods,str):
        methods = [methods]

    dsr = ds.to_dataframe().resample(freq)
    dsouts = []
    for method in methods:
        # what we want (quickly), but in Pandas form
        df_h = dsr.apply(method)
        # rebuild xarray dataset with attributes
        vals = []
        for c in df_h.columns:
            vals.append(
                xr.DataArray(data=df_h[c],
                             dims=['time'],
                             coords={'time': df_h.index},
                             attrs=ds[c].attrs)
            )
        dsouts.append(xr.Dataset(dict(zip(df_h.columns, vals)), attrs=ds.attrs))

    if len(dsouts) == 1:
        dsouts = dsouts[0]
    return dsouts

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

def get_sensor_config(sconfig: dict|None = None) -> dict:
    """ Read the sensor configuration from the default json file and merge if needed. 
    """
    fn_config = pkg_res.resource_filename("pyrnet", "share/pyrnet_sensor_config.json")
    default_config = pyrutils.read_json(fn_config)
    if sconfig is None:
        sconfig = default_config
    sconfig = {**default_config, **sconfig}
    return sconfig
    
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

def calc_encoding(sconfig:dict, ADCV=3.3, ADCbits=10) -> dict:
    ADCfac = ADCV / (2**ADCbits-1) # Last bit is reserved 
    sencoding ={}
    for k,v in sconfig.items():
        sencoding.update(
            {k: dict(
                units=v['units'],
                scale_factor=v['C']*ADCfac/v['gain'],
                add_offset=v['offset'],
                valid_range= [0, min(((2**ADCbits-1), int(v['Vmax']*v['gain']/ADCfac)))]
            )}
        )
    return sencoding


# %% ../../nbs/pyrnet/data.ipynb 18
def add_encoding(ds, vencode=None):
    """
    Set valid_range attribute and encoding to every variable of the dataset.

    Parameters
    ----------
    ds: xr.Dataset
        Dataset of any processing level. The processing level will be
        determined by the global attribute 'processing_level'.
    vencode: dict or None
        Dictionary of encoding attributes by variable name, will be merged with pyrnet default cfmeta. The default is None.

    Returns
    -------
    xr.Dataset
        The input dataset but with encoding and valid_range attribute.
    """
    # prepare netcdf encoding
    _, vattrs_default, vencode_default = get_cfmeta()

    # Add valid range temporary to encoding dict.
    # As valid_range is not implemented in xarray encoding,
    #  it has to be stored as a variable attribute later.
    for k in vencode_default:
        if "valid_range" not in vencode_default[k]:
            continue
        vencode_default = assoc_in(vencode_default,
                                   [k,'valid_range'],
                                   vattrs_default['valid_range'])

    # merge input and default with preference on input
    if vencode is None:
        vencode = vencode_default
    else:
        a = vencode_default.copy()
        b = vencode.copy()
        vencode = {}
        for k in set(a)-set(b):
            vencode.update({k:a[k]})
        for k in set(a)&set(b):
            vencode.update({k: {**a[k],**b[k]}})
        for k in set(b)-set(a):
            vencode.update({k:b[k]})

    # add encoding to Dataset
    for k, v in vencode.items():
        if k not in ds.keys():
            continue
        for ki in [key for key in ds if key.startswith(k)]:
            ds[ki].encoding = v
        if "valid_range" not in vencode[k]:
            continue
        # add valid_range to variable attributes
        for ki in [key for key in ds if key.startswith(k)]:
            ds[ki].attrs.update({
                'valid_range': vencode[k]['valid_range']
            })

    # special treatment of time and flux variables
    if ds.processing_level=='l1a':
        ds["gpstime"].encoding.update({
            "dtype": 'f8',
            "units": f"seconds since {np.datetime_as_string(ds.gpstime.data[0], unit='D')}T00:00Z",
        })
        ds["adctime"].encoding.update({
            "units": f"milliseconds",
        })
    elif ds.processing_level == 'l1b':
        # ds = stretch_resolution(ds)
        # # special treatment for flux variables
        # for k in ['ghi', 'gti']:
        #     if k not in ds:
        #         continue
        #     # add encoding
        #     dtype = ds[k].encoding['dtype']
        #     int_limit = np.iinfo(dtype).max
        #     valid_range = [0, int_limit - 1]
        #     scale_factor = 1500. / float(int_limit - 1)
        #     ds[k].encoding.update({
        #         "scale_factor": scale_factor,
        #         "_FillValue": int_limit,
        #     })
        #     ds[k].attrs.update({
        #         "units": "W m-2",
        #         "valid_range": valid_range
        #     })
        ds["time"].encoding.update({
            "dtype": 'f8',
            "units": f"seconds since {np.datetime_as_string(ds.time.data[0], unit='D')}T00:00Z",
        })
    else:
        raise ValueError("Dataset has no 'processing_level' attribute.")
    return ds

# %% ../../nbs/pyrnet/data.ipynb 22
def to_l1a(
        fname : str,
        *,
        station: int,
        report: dict|pd.DataFrame|None,
        date_of_measure : np.datetime64 = np.datetime64("now"),
        config: dict|None = None,
        sconfig: dict|None = None,
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
    sconfig: dict
        Config for ADC and amplifier for each sensor. The default is "../share/pyrnet_sensor_config.json"
    global_attrs: dict
        Additional global attributes for the Dataset. (Overrides cfmeta.json attributes)
    Returns
    -------
    xarray.Dataset
        Raw Logger data for one measurement periode.
    """
    ADCV = 3.3
    ADCbits = 10
    
    # load and merge  default config
    config = get_config(config)
    gattrs, vattrs, vencode = get_cfmeta(config)
    
    # update encoding with sensor config for l1a
    sconfig = get_sensor_config(sconfig)
    sencoding = calc_encoding(sconfig, ADCV=ADCV, ADCbits=ADCbits)
    # update encoding with sensor config for l1a
    for var, enc in sencoding.items():
        for k, v in enc.items():
            if k in ["valid_range", "units"]:
                vattrs = assoc_in(vattrs, [var, k], v)
            else:
                vencode = assoc_in(vencode, [var, k], v)
                
    # update additional global attributes
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
    adc_volts = ADCV * rec_adc[:,2:] / float(2**ADCbits - 1)

    # 2. Get Logbook maintenance quality flags
    key = f"{station:03d}"
    if report is None:
        logger.warning("No report available!")
        report = {}
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
    vattrs = assoc_in(vattrs, ["gti","hangle"], 0.)
    vattrs = assoc_in(vattrs, ["gti","vangle"], 0.)
    # update with angles from mapping file
    if config['gti_angles'] is not None:
        gti_angles = pyrutils.read_json(config['file_gti_angles'])[config['gti_angles']]
        if key in gti_angles:
            hangle = np.nan if gti_angles[key][0] is None else gti_angles[key][0]
            vangle = np.nan if gti_angles[key][1] is None else gti_angles[key][1]
            vattrs = assoc_in(vattrs, ["gti","hangle"], hangle)
            vattrs = assoc_in(vattrs, ["gti","vangle"], vangle)

    if adc_volts.shape[1]<5: # gti data is not available
        adc_volts = np.concatenate((adc_volts,-1*np.ones(adc_volts.shape[0])[:,None]),axis=1)

    # 8. Make xarray Dataset
    values = {}
    for k,v in sconfig.items():
        offset = sconfig[k]["offset"]
        gain = sconfig[k]["gain"]
        C = sconfig[k]["C"]
        iadc = sconfig[k]["iadc"]
        volts = adc_volts[:,iadc][:,None]
        values.update(
            {k: offset + C*volts/gain}
        )
    ds = xr.Dataset(
        data_vars={
            "ghi": (("adctime","station"), values["ghi"]), # [V]
            "gti": (("adctime","station"), values["gti"]), # [V]
            "ta": (("adctime","station"), values["ta"]), # [K]
            "rh": (("adctime","station"), values["rh"]), # [-]
            "battery_voltage": (("adctime","station"), values["battery_voltage"]), # [V]
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
    ds = add_encoding(ds, vencode)

    return ds

# %% ../../nbs/pyrnet/data.ipynb 55
def to_l1b(
        fname: str,
        *,
        config: dict | None = None,
        global_attrs: dict | None = None,
        check_adc_sync: bool = True
) -> xr.Dataset|None:

    config = get_config(config)
    gattrs, vattrs, vencode = get_cfmeta(config)

    if global_attrs is not None:
        gattrs.update(global_attrs)

    # 1. Load l1a data
    ds_l1a = xr.open_dataset(fname)
    # check correct file
    if ds_l1a.processing_level != "l1a":
        logger.warning(f"{fname} is not a l1a file. Skip.")
        return None

    # 2. Sync GPS to ADC time
    adctime = pyrlogger.sync_adc_time(
        adctime = ds_l1a.adctime.values,
        gpstime = ds_l1a.gpstime.values,
        iadc = ds_l1a.iadc.squeeze().values.astype(int),
        check_results = check_adc_sync
    )
    
    if adctime is None:
        logger.warning(f"Could not fit GPS to ADC time for file {fname}. Skip.")
        return None

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
    if (ds_l1b.time.values[0] + 3*stripminutes) > ds_l1b.time.values[-1]:
        logger.warning(f"{fname} has not enough data. Skip.")
        return None

    ds_l1b = ds_l1b.isel(time=ds_l1b.time>ds_l1b.time.values[0] + stripminutes)
    ds_l1b = ds_l1b.isel(time=ds_l1b.time<ds_l1b.time.values[-1] - stripminutes)
    if ds_l1b.time.size < 10:
        logger.warning(f"{fname} has not enough data, after strip. Skip.")
        return None


    logger.info(f"Dataset time coverage after strip: {ds_l1b.time.values[0]} - {ds_l1b.time.values[-1]}")

    # 5. rad flux calibration
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


    # calibrate radflux variables
    for i, radflx in enumerate(config['radflux_varname']):
        if cfac[i] is None:
            # drop if calibration/instrument don't exist (probably secondary pyranometer).
            ds_l1b = ds_l1b.drop_vars([var for var in ds_l1b if radflx in var])
            continue
        ds_l1b[radflx].values = ds_l1b[radflx].values*1e6/(cfac[i]) # V -> W m-2
        ds_l1b[radflx].attrs.update({
            "serial": serial[i],
            "calibration_factor": cfac[i],
            "calibration_function": "flux (W m-2) = flux (V) * calibration_factor (W m-2 V-1)",
        })


    # 6. resample to desired resolution
    # save station coordinate
    station_dim = {"station": ds_l1b["station"].values}
    
    # resample on time dimension with specified methods
    methods = ['mean'] + config["l1b_resample_stats"]
    res = resample(
        ds_l1b.squeeze().drop_vars("station"), # drop station coordinate and variable
        freq=config['l1bfreq'],
        methods=methods,
        kwargs=dict(skipna=True)
    )
    
    # add standard names for new variables
    ds_l1b = res[0]
    for i, method in enumerate(methods[1:]):
        for var in config["radflux_varname"]:
            ds_l1b[f"{var}_{method}"] = res[i+1][var]
            ds_l1b[f"{var}_{method}"].attrs.update({
                "standard_name": f"{method}_"+ds_l1b[f"{var}_{method}"].attrs["standard_name"]
            })
    
    # add station dimension back again
    ds_l1b = ds_l1b.expand_dims(station_dim, axis=-1)
    
    # stretch valid range to not lose resolution due to averaging
    # TODO: remove stretch resolution 
    # ds_l1b = stretch_resolution(ds_l1b)

    # 7. Interpolate GPS coordinates to l1b time
    ds_gps = ds_l1a.drop_dims("adctime")
    ds_gps = ds_gps.drop_vars(['iadc'])

    # Decide whether geo coordinates should be averaged or not
    if config['average_latlon']:
        ds_gps = ds_gps.mean('gpstime', skipna=True, keep_attrs=True)
    else:
        ds_gps = ds_gps.interp(gpstime=ds_l1b.time)
        ds_gps = ds_gps.drop_vars("gpstime")

    ds_l1b = xr.merge((ds_l1b,ds_gps))

    # 8. Calc and add sun position
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

    # add global coverage attributes
    ds_l1b = update_coverage_meta(ds_l1b, timevar="time")
    ds_l1b.attrs["processing_level"] = 'l1b'
    now = pd.to_datetime(np.datetime64("now"))
    ds_l1b.attrs["history"] = ds_l1b.history + f"{now.isoformat()}: Generated level l1b  by pyrnet version {pyrnet_version}; "

    # update encoding
    ds_l1b = add_encoding(ds_l1b, vencode=vencode)

    return ds_l1b
