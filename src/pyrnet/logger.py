# AUTOGENERATED! DO NOT EDIT! File to edit: ../../nbs/pyrnet/logger.ipynb.

# %% auto 0
__all__ = ['logger', 'dtype_gprmc', 'parse_gprmc', 'parse_adc', 'read_records', 'get_adc_time', 'sync_adc_time', 'adc_binning',
           'resample_mean', 'interpolate_coords']

# %% ../../nbs/pyrnet/logger.ipynb 2
from numpy.typing import NDArray,ArrayLike
import re
import gzip
import numpy as np
import pandas as pd
from scipy.stats import linregress
import logging

from . import utils

# logging setup
logging.basicConfig(
    filename='pyrnet.log',
    encoding='utf-8',
    level=logging.DEBUG,
    format='%(asctime)s %(name)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)


# Local variables
_nat      = np.datetime64('NaT','ms')
_re_gprmc = re.compile('^(\d+,\d+ \d) \$GPRMC,(.*)(\*\w{2})$')
_re_adc   = re.compile('^(\d+)(\s\d+)+$')

# %% ../../nbs/pyrnet/logger.ipynb 11
def parse_gprmc(s, date_of_measure=np.datetime64('now')):
    '''
    Parse a string with a GPRMC GPS record

    Parameters
    ----------
    s: string
        A GPRMC record
    date_of_measure: datetime or datetime64
        A rough time, when the measurements happen to account for GPS rollover.
        Precise datetime is only necessary for the period of 2019-05-06 to 2019-08-17.
        Otherwise, providing a year is sufficient.
        The default is np.datetime64("now"), which is sufficient for all measurements conducted later than 2019-08-17.

    Returns
    -------
    gprmc : tuple
        A tuple with datetime64,status,lat,lon
    '''
    date_of_measure = utils.to_datetime64(date_of_measure)
    # split fields of GPRMC record
    f = s.split(',')
    status = f[1]
    if status=='A':
        try:
            # parse latitude
            lat = float(f[2][:2])+float(f[2][2:])/60
            if f[3]=='S':
                lat *= -1.0
            # parse longitude
            lon = float(f[4][:3])+float(f[4][3:])/60
            if f[5]=='W':
                lon *= -1.0
            if date_of_measure>np.datetime64("2019-04-06") and date_of_measure<np.datetime64("2019-08-17"): #account for gps week rollover
                YY = '19'+f[8][4:6]
            else:
                YY = "20"+f[8][4:6]
            mm = f[8][2:4]
            dd = f[8][0:2]
            HH = f[0][0:2]
            MM = f[0][2:4]
            SS = f[0][4:]
            dt = np.datetime64(YY+'-'+mm+'-'+dd+'T'+HH+':'+MM+':'+SS,'ms')
            if date_of_measure>np.datetime64("2019-04-06"): #account for gps week rollover -> date jump 1024 weeks back at 2019-04-06
                dt= dt+np.timedelta64(1024,'W')
            r = (dt,status,lat,lon)
        except:
            return (_nat,'V',np.nan,np.nan)
    else:
        r = (_nat,'V',np.nan,np.nan)
    return r

# %% ../../nbs/pyrnet/logger.ipynb 13
def parse_adc(s):
    '''
    Parse an ADC record

    Parameters
    ----------
    s: string
        The ADC record

    Returns
    -------
    t: tuple
       A tuple of digital counts of the ADC
    '''
    return tuple(map(int,s.split()))

# %% ../../nbs/pyrnet/logger.ipynb 15
dtype_gprmc = [
    ( 'time',   'datetime64[ms]' ),
    ( 'status', 'S1' ),
    ( 'lat',    'f8' ),
    ( 'lon',    'f8' ),
    ( 'iadc',   'u4' )
]

def read_records(fname: str,
                 date_of_measure: np.datetime64 = np.datetime64('now')) -> (NDArray, NDArray):
    '''
    Read the GPRMC and ADC records from the pyranometer logger files

    Parameters
    ----------
    fname: string
        The filename of the logger file
    date_of_measure: numpy.datetime64
        Date of measurement to account for gps rollover

    Returns
    -------
    rec_adc: ndarray
        The 10bit ADC readings
    rec_gprmc: recarray
        The GPRMC GPS records
    '''
    logger.info(f"Start reading records from file: {fname}")
    date_of_measure = utils.to_datetime64(date_of_measure)
    # Read file, use errors='ignore' to skip non UTF-8 characters
    # non UTF-8 characters may arise in broken GPS strings from time to time
    if fname[-3:]=='.gz':
        f = gzip.open(fname,'rt',errors='ignore') # open in text mode, ignore non UTF8 characters
    else:
        f = open(fname,'r',errors='ignore')
    lines =[l.rstrip() for l in f.readlines()]
    f.close()

    ##- skip almost empty files
    if len(lines)<20:
        logger.info("Skip file, as number of records is < 20.")
        return False,False

    # remove last line -> mostly damaged or empty
    lines=lines[:-1]
    # remove gps line at the end -> else processing issues
    if _re_gprmc.match(lines[-1]):
        lines=lines[:-1]

    rec_gprmc = []
    rec_adc = []
    iadc = 0
    for i,l in enumerate(lines):
        m = _re_gprmc.match(l)
        if m:
            r = parse_gprmc(m.group(2), date_of_measure)
            if not np.isnat(r[0]):
                # add number of adc values before GPS line
                rec_gprmc.append(r+(iadc,))
        elif _re_adc.match(l):
            r = parse_adc(l)
            if iadc==0:
                adc_len=len(r)
            # if record line is incomplete (due to power cut off)
            # the line is dropped
            if len(r)==adc_len:
                rec_adc.append(r)
                iadc += 1
        else:
            # unhandled record...
            pass
    rec_adc   = np.array(rec_adc,dtype=np.uint16)
    rec_gprmc = np.array(rec_gprmc,dtype=dtype_gprmc).view(np.recarray)
    logger.info("Done reading records from raw file.")
    return rec_adc, rec_gprmc

# %% ../../nbs/pyrnet/logger.ipynb 20
def get_adc_time(rec_adc):
    """
    Get Milliseconds from Start of ADC measurement.

    Parameters
    ----------
    rec_adc: ndarray
        The digital counts of the ADC from the logger file.

    Returns
    -------
    ndarray
        Milliseconds from start of the measurement.
    """
    # get millisecond part
    ta = rec_adc[:,0].astype('timedelta64[ms]')
    # get time difference between records
    dt = np.diff(ta)
    dt[dt<np.timedelta64(-850,'ms')] += 1000
    # get cummulative time offset rel. to first ADC record
    ta[0] = 0
    ta[1:] = np.cumsum(dt)
    return ta

# %% ../../nbs/pyrnet/logger.ipynb 24
def sync_adc_time(adctime, gpstime, iadc, check_results=True):
    '''
    Synchronize the ADC time to the GPS records

    Parameters
    ----------
    adctime: ndarray
        Milliseconds from start of the ADC measurement.
    gpstime: ndarray
        GPS time
    iadc: ndarray of int
        Index of the last ADC sample before a GPS record has been stored.
    check_results: bool
        If True, check plausibility of fitted slope and offset (abs(slope)<10s/day; abs(offset)<2s).
        The default is True

    Returns
    -------
    time: ndarray(datetime64[ms])
        The time of the ADC records
    '''
    # assure milliseconds
    ta = adctime.astype('timedelta64[ms]')
    # assure int type
    iadc = np.array(iadc).astype(int)
    # get GPStime of the previous ADC record
    # here we assume, that the last ADC sample stored before a GPS record
    # has the same time as the GPS timestamp
    t1 = ta[iadc]/np.timedelta64(1,'ms')
    # time of GPS records from GPRMC record
    t2 = (gpstime-gpstime[0])/np.timedelta64(1,'ms')
    a,b = linregress(t1,t2)[:2]
    t = gpstime[0]+ta*a+b.astype('timedelta64[ms]')

    drift = (1/a-1)*86400
    logger.info('Sync ADC time to GPS Fit Summary:')
    logger.info('|-- Drift  : {0:7.2f} [s/day]'.format(drift))
    logger.info('|-- Slope  : {0:13.8f}'.format(a))
    logger.info('|-- Offset : {0:7.2f} [s]'.format(b/1000))
    logger.info('|-- Jitter : {0:7.2f} [ms]'.format(np.std(t2-(a*t1+b))))
    
    if check_results:
        if np.abs(drift)>10:
            logger.warning("Absolute ADC drift larger than 10 s/day.")
            return None
    return t

# %% ../../nbs/pyrnet/logger.ipynb 36
def adc_binning(rec_adc, time, bins=86400):
    """
    Binning and averaging of ADC samples

    Parameters
    ----------
    rec_adc: ndarray
        The ADC records parsed from the logger.
    time: ndarray of time objects
        Sample time of ADC records.
    bins: int
        Number of desired bins per day. The default is 86400, which result in
        mean values of 1 second steps per day. Maximum resolution is 86400000.

    Returns
    -------
    ndarray, ndarray(datetime64)
        Binned ADC records and corresponding time.
    """
    time = utils.to_datetime64(time)
    # starting day
    t0 = time[0].astype('datetime64[D]')
    # convert time to 'days from t0'
    dday = (time-t0)/np.timedelta64(1,'D')
    # calculate time bins of output dataset
    it = np.int64(dday*bins)
    # index for unique bins (inv_idx) and count of samples per bin (cnt)
    uval, inv_idx, cnt = np.unique(it,
                                   return_inverse=True,
                                   return_counts=True)
    logger.info(f"ADC records fill {len(uval)} bins of data.")
    # Calculate average of sample values per bin
    # The first two columns of rec_adc will be omitted as they store the
    # internal measures for timing and battery (first two columns)
    V = np.zeros((len(uval),rec_adc.shape[1]-2))
    for i in range(V.shape[1]):
        V[:,i] = np.bincount(inv_idx,weights=rec_adc[:,i+2])/cnt
    bintime = t0+ np.timedelta64(86400000,'ms')*uval.astype(np.float64)/bins
    logger.info(f"ADC records span a time period from {bintime[0]} to {bintime[-1]}.")
    return V, bintime

# %% ../../nbs/pyrnet/logger.ipynb 38
def resample_mean(ds,freq='1s'):

    # start and end bin time
    start_time = np.datetime64(
       pd.to_datetime(ds.time.values[0]).floor(freq)
    )
    end_time = np.datetime64(
        pd.to_datetime(ds.time.values[-1]).floor(freq)
    )

    # bin time
    bintime = pd.date_range(
        start_time,
        end_time,
        freq=freq
    ).floor(freq)

    ds_r = ds.assign_coords(
        {
            "time_resampled": ("time_resampled", bintime)
        }
    )

    # calculate bin index of output dataset
    it = np.int64(
        (ds.time.values - start_time)/pd.Timedelta(freq)
    )

    # index for unique bins (inv_idx) and count of samples per bin (cnt)
    uval, inv_idx, cnt = np.unique(it,
                                   return_inverse=True,
                                   return_counts=True)

    # apply to all time dependent variables
    for var in ds:
        if 'time' in ds[var].dims:
            # replace time dimension with time_resampled
            vardims = ds[var].dims
            newdims = [d if d!='time' else 'time_resampled' for d in vardims]
            # if only time dimension, take the shortcut
            if len(vardims)==1:
                newval = np.bincount(inv_idx,weights=ds[var].values)/cnt
                ds_r = ds_r.assign( {var: (newdims, newval)})
            elif len(vardims)==2: # more than the time dimension, assume 2
                N = ds[var].shape[1]
                newval = np.zeros((bintime.size,N))
                for i in range(N):
                    newval[:,i] = np.bincount(inv_idx,weights=ds[var].values[:,i])/cnt
                ds_r = ds_r.assign( {var: (newdims, newval)})
            else:
                raise ValueError("logger.resample is implemented for 2dims only.")
        # add attributes again
        ds_r[var].attrs.update(ds[var].attrs)
        ds_r[var].encoding.update(ds[var].encoding)

    # drop original time and rename
    ds_r = ds_r.drop_dims("time").rename({"time_resampled":"time"})
    # add time encoding
    ds_r.time.encoding.update({
        "units": f"seconds since {np.datetime_as_string(ds_r.time.data[0], unit='D')}T00:00Z",
    })

    return ds_r


# %% ../../nbs/pyrnet/logger.ipynb 40
def interpolate_coords(rec_gprmc, time):
    """
    Interpolate lat and lon from gps records

    Parameters
    ----------
    rec_gprmc: recarray
        The GPRMC records from the logger file
    time: list or array of time objects

    Returns
    -------
    ndarray, ndarray
        lat and lon interpolated from GPS records to `time`
    """
    time = utils.to_datetime64(time)
    t0 = time[0].astype('datetime64[D]')
    # coordinate variables
    x1 = (time-t0)/np.timedelta64(1,'ms')
    x2 = (rec_gprmc.time-t0)/np.timedelta64(1,'ms')
    lat = np.interp(x1,x2,rec_gprmc.lat)
    lon = np.interp(x1,x2,rec_gprmc.lon)
    return lat, lon
