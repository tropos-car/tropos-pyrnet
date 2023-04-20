import pyproj
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

# python -m pip install git+https://github.com/hdeneke/trosat-base.git#egg=trosat-base
from trosat import sunpos as sp



campaign_pfx = {
    'eifel': 'hope',
    'hope_juelich': 'hope',
    'hope_melpitz': 'hopm',
    'lindenberg': 'ioprao',
    'melcol': 'mcol',
}
    
data_dir = "https://tds.tropos.de/thredds/dodsC/scccJher/{dt:%Y}_{campaign}/old/nc/"

fname_fmt = '{campaign_pfx}_trop_pyrnet00_l1_rsds_v00_{dt:%Y%m%d}000000.nc'



solconst = 1359.0
max_missing = 300       # Maximum allowed number of missing records
min_good    = 85400     # Minimum allowed number of good records

geod       = pyproj.Geod(ellps='WGS84')

def read_hdcp2( dt, fill_gaps=True, resample=False, campaign='hope_juelich' ):
    '''
    Read HDCP2-formatted datafiles from the pyranometer network
    
    Arguments
    ---------
    dt: datetime.date
        The date of the data to read
    fill_gaps: bool
        A flag indicating whether gaps should be filled by interpolation
    resample: bool
        not implemented, no effects
    campaign: str
        specify campaign ['eifel','hope_juelich','hope_melpitz','lindenberg','melcol']

    Returns
    -------
    dataset : xarray.Dataset 
        The pyranometer network observations
    '''
    # load dataset
    fname = data_dir+fname_fmt
    fname = fname.format(dt=dt,
                         campaign=campaign,
                         campaign_pfx=campaign_pfx[campaign])
    ds = xr.open_dataset(fname,mask_and_scale=False)

    # select good stations
    igood = (np.sum(ds.rsds.data<-900.0,axis=0)<1000)&(np.sum(ds.rsds_flag==1,axis=0)>85400)
    ds = ds.isel(nstations=igood)

    # fill gaps if requested
    if fill_gaps==True:
        x = (ds.time-ds.time[0])/np.timedelta64(1,'s')
        for i in np.arange(ds.dims['nstations']):
            y = ds.rsds.data[:,i]
            m = y>-990.0
            if not np.all(m):
                f = interp1d(x[m],y[m],'linear',bounds_error=False,fill_value='extrapolate')
                ds.rsds[~m,i]=f(x[~m])
    # add additional DataArrays
    jd = (ds.time.data-np.datetime64(sp.EPOCH_JD2000_0))/np.timedelta64(1,'D')
    ds['esd'] = sp.earth_sun_distance(jd[0]+0.5)
    szen = sp.sun_angles(jd[:,None],ds.lat.data[None,:],ds.lon.data[None,:])[0]
    ds['szen']    = xr.DataArray(szen,dims=('time','nstations'),coords={'time':ds.time.data})
    ds['mu0']     = np.cos(np.deg2rad(ds.szen))
    ds['gtrans']  = ds.rsds/ds.esd**2/solconst/ds['mu0']
    return ds.rename({'rsds_flag':'qaflag','rsds':'ghi'})    

def get_xy_coords(lon, lat, lonC=None, latC=None):
    '''
    Calculate Cartesian coordinates of network stations, relative to the mean
    lon/lat of the stations
    '''
    n  = len(lon)
    if not lonC:
        lonC = lon.mean()
    if not latC:
        latC = lat.mean()
    x = np.zeros(n)
    y = np.zeros(n)
    az,_,d = np.array([geod.inv(lonC, latC, lon[i], lat[i]) for i in range(n)]).T
    x = d*np.sin(np.deg2rad(az))
    y = d*np.cos(np.deg2rad(az))
    return x,y

def pairwise_distance_matrix( x, y ):
    '''
    Get square matrix with Euclidian distances of stations
    '''
    return np.sqrt( (x[None,:]-x[:,None])**2+(y[None,:]-y[:,None])**2 )
