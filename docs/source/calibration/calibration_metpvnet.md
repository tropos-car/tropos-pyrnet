
# 2019-06 MetPVNet
We calibrate the Pyrnet (30 stations with two pyranometer setup) for MetPVNet 2018/2019 measurement campaigns.
Cross-calibration is done versus MORDOR3 global pyranometer.

Pyranometer data is processed with Pyrnet processing software: https://gitea.tropos.de/walther/Pyrnet_processing.git (20.06.2019)


```python
#[notes]
#Python cells begin with a tag and dependency line: 
#[tag] =>[tag_of_cell_mandatory_to_run_before]
```

## Imports 


```python
#[imports]
# python modules
import os
import xarray as xr
import numpy as np
import pandas as pd
from matplotlib import pyplot
```

## Definitions
Set data paths and threshholds


```python
#[definitions] =>[imports]
# database path 
# mordor path suffix - with wildcards to get all files for that period
mpf="mordor/201906_tropos/ncdata/2019/06/*_Rad*.nc"

# pyrnet path suffix - lvl1 means uncalibrated raw data in nc standard format
ppf="/pyrnet/201906_troposcalib/lvl1/"

#get involved station IDs of pyranometer net
stations=[]
for fname in os.listdir(pf+ppf):
    if fname[:3]=='Pyr':
        stations.append(int(fname.split('_')[1]))
stations=np.sort(stations)
#stations=[1]
# Factor for ADC gain (300) and unit transformation [V -> uV]
gain_units_factor=1e6/300.0 
        
#calibration periode
start_date=np.datetime64("2019-06-08")
end_date=np.datetime64("2019-06-17")
```

## MORDOR3 dataset
Currently (07.2019) MORDOR netcdf data is still not fully developed in terms of variable name / dimension connection for optimal xarray use. This will be change in future. Check the code carefully if you run for later calibration


```python
#[mordor_data] =>[definitions]
#load mordor radiation data
mds=xr.open_mfdataset(pf+mpf)
#adjust handling of dimension and time variable in mordor nc file
mds=mds.rename({'t':'time'})
mds=mds.set_coords('time')
## Preselect data to reduce ram usage
#drop data after pyrnet was removed from roof
mds=mds.where(mds['time']<end_date,drop=True)
#drop data of first two days because of time synchronous issues of measure laptop
mds=mds.where(mds['time']>start_date,drop=True)

#print mds header
print(mds)




```

    <xarray.Dataset>
    Dimensions:  (time: 7749295)
    Coordinates:
      * time     (time) datetime64[ns] 2019-06-08T00:00:00.100000 ... 2019-06-16T23:59:59.900000
    Data variables:
        GHI      (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
        UPI_lw   (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
        UPI      (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
        DNI      (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
        szen     (time) float32 dask.array<shape=(7749295,), chunksize=(861142,)>
        sazi     (time) float32 dask.array<shape=(7749295,), chunksize=(861142,)>
        DHI      (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
        GHI_lw   (time) float64 dask.array<shape=(7749295,), chunksize=(861142,)>
    Attributes:
        earth_sun_distance:        1.01490394190788
        earth_sun_distance_units:  AU
        Title:                     Radiaton data observed with MORDOR3
        Institution:               Leibniz Institute for Tropospheric Research (T...
        Contact_person:            Hartwig Deneke plus satellite group, sat@tropo...
        Source:                    Product derived from MObile RaDiation ObseRvat...
        History:                   Data file generated at TROPOS_calibration at 2...
        Conventions:               MESOR
        Author:                    Jonas Witthuhn (witthuhn@tropos.de)
        License:                   For non-commercial use only.
        latitude:                  51.35
        longitude:                 12.44
        latitude_units:            deg N
        lontitude_units:           deg E


## Calibration
Calibration follows ISO 9847 in <cite>[Indian Standard - Calibration Field Pyranometers by comparison to a reference pyranometer] [1]</cite>.

Cloudy sky threatment is applied due to variing conditions during the calibration periode.

### Step 1 
Drop Night measures and low signal measures from pyranometer data. Since calibration without incoming radiation doesnt work.

These data is kept for calibration:

* solar zenith angle < 80° ( as recommended in <cite>[1]</cite>)
* Measured Voltage > 0.033 V (drop the lowest 1%)

### Step 2

Fit MORDOR reference data to selected Pyranometer data.

### Calculating the calibration factor
We define here the equations used before moving to step 3.

Calculate series of calibration factor for this dataset following <cite>[1]</cite> 5.4.1.1 equation (1).

Our final calibration equation looks like this:
$$ Irradiance = Volts * 10^{-6} * gain^{-1} * calibration^{-1} $$
With:
* Irradiance: Calibrated irradiance Value [Wm-2]
* Volts: Measured voltage [V] -> to [uV] with $10^{-6}$
* gain: ADC gain = 300
* calibration: calibration factor with units [uV W-1 m2]
Therefore we calculate:
$$ C(ij) = \frac{V(ij) * 10^{-6}}{300 * I(ij)} $$
, for each timestep (j, hourly series ; i instant sample).
* I - Irradiance from MORDOR
* V - Voltage from Pyranometer to calibrate
* C - resulting calibration factor

To get rid of outliers/ arbitrary shading etc. we calculate the integratet series of calibration following <cite>[1]</cite> 5.4.1.1 equation (2):
$$ C(j) = \frac{10^{-6}}{300}*\frac{[V(j)]_{int}}{[I(j)_{int}]} $$
, with $ []_{int} $ beeing integrated values over intervals of one hour.

### Step 3
Remove outliers from series using xarray grouping and apply function. The following functions removes outliers (deviation more than 2% according to <cite>[1]</cite>) from a selected group. This step involves calculating calibration series and the integration of one hour intervals to smooth out high variable situation, which would break the calibration even when time synchronization is slightly off. Also this gets rid of some random shading events like birds / chimney / rods in line of sigth, which would affect calibration otherwise. We following <cite>[1]</cite> 5.4.1.1 equation (2) here. 

### Step 4
The series of measured voltage and irradiance is now without outliers. So we use equation 1 again to calculate from this reduced series the calibration factor for the instant samples.

### Step 5
 We just found the Calibration factor to be the mean of the reduced calibration factor series and the uncertainty to be the standard deviation of this reduced series. Steo 3, 4 and 5 are done for every pyranometer seperate

    
    [1] https://archive.org/details/gov.in.is.iso.9847.1992


```python
#[remove_function] => [imports]
##function for outlier removal
def remove_outliers(x):
    """
    x is an xarray dataset containing these variables:
    coords: 'time' - datetime64
    'gain_units_factor' - float - to scale the calibration
    'volts' - array - voltage measures of pyranometer
    'irradiance' - array - measured irradiance of mordor
    """
    #calculate calibration series for single samples
    C=x['volts']*x['gain_units_factor']/x['irradiance']
    #integrated series 
    ix=x.integrate('time')
    M=ix['volts']*x['gain_units_factor']/ix['irradiance']
    while np.any(np.abs(C-M)>0.02*M):
        #calculate as long there are outliers deviating more than 2 percent
        x=x.where(np.abs(C-M)<0.02*M)
        C=x['volts']*x['gain_units_factor']/x['irradiance']
        #integrated series 
        ix=x.integrate('time')
        M=ix['volts']*x['gain_units_factor']/ix['irradiance']
    #return the reduced dataset x
    return x
        
```


```python
#[calibrations] =>[mordor_data][remove_function]

#preparing output table
print('BOX , C(pyr1), stdC(pyr1) , C(pyr_tilt) , stdC(pyr_tilt)')
with open('new_calibration.csv','w') as out:
    out.write('#BOX , C(pyr1), stdC(pyr1) , C(pyr_tilt) , stdC(pyr_tilt)\n')
    
# do calibration for all stations
for st in stations:
    #load uncalibrated pyrnet data
    pds=xr.open_mfdataset(str(pf+ppf+
                              "Pyr_%s/"%(str(st).zfill(3))+
                              "metpvnet_preliminary_PYRNET_*_%s_lvl1.nc"%(str(st).zfill(3))))
    #drop some nan values along time dimension
    pds=pds.dropna('time')
    
    ## Step 1 #############################################
    #drop data after pyrnet was removed from roof
    pds=pds.where(pds['time']<end_date,drop=True)
    #drop data of first two days because of time synchronous issues of measure laptop
    pds=pds.where(pds['time']>start_date,drop=True)
    #drop night measurements
    pds=pds.where(pds['sza']<80,drop=True)
    #drop low and nan measures
    pds=pds.where(pds['ghi']>0.033,drop=True)
    
    ## Step 2 #############################################
    #fit mordor measures to pyranometer
    tds=mds.reindex_like(pds)
    
    ## Step 3 #############################################
    ##Calibration pyr1
    #preparing selection dataset for calibration
    Cds=xr.Dataset({'time':pds['time'],
                    'gain_units_factor':gain_units_factor,
                    'irradiance':('time',tds['GHI']),
                    'volts':('time',pds['ghi'])})
    #do hourly samples  as reccomended in ISO 9842
    Cds=Cds.groupby('time.hour').apply(remove_outliers)
    
    ## Step 4 #############################################
    # calculate the calbration series without outliers
    C=Cds['volts']*gain_units_factor/Cds['irradiance']
    
    ## Step 5 #############################################
    # we found the calibration for pyranometer 1
    Cpyr1=np.nanmean(C.values)
    uCpyr1=np.nanstd(C.values)
    
    
    ## Step 3 #############################################
    ##Calibration pyr2
    #preparing selection dataset for calibration
    Cds=xr.Dataset({'time':pds['time'],
                    'gain_units_factor':gain_units_factor,
                    'irradiance':('time',tds['GHI']),
                    'volts':('time',pds['ghi_tilt'])})
    #do hourly samples  as reccomended in ISO 9842
    Cds=Cds.groupby('time.hour').apply(remove_outliers)
    
    ## Step 4 #############################################
    # calculate the calbration series without outliers
    C=Cds['volts']*gain_units_factor/Cds['irradiance']
    
    # Step 5 #############################################
    # we found the calibration for pyranometer 2
    Cpyr2=np.nanmean(C.values)
    uCpyr2=np.nanstd(C.values)
    
    print("%d , %.2f , %.3f , %.2f , %.3f"%(st,Cpyr1,uCpyr1,Cpyr2,uCpyr2))
    with open('new_calibration.csv','a') as out:
        out.write("%d , %.2f , %.3f , %.2f , %.3f\n"%(st,Cpyr1,uCpyr1,Cpyr2,uCpyr2))
```

    BOX , C(pyr1), stdC(pyr1) , C(pyr_tilt) , stdC(pyr_tilt)
    1 , 7.73 , 0.244 , 6.98 , 0.211
    4 , 7.62 , 0.462 , 8.23 , 0.247
    5 , 7.48 , 0.259 , 6.37 , 0.224
    7 , 7.52 , 0.242 , 7.20 , 0.277
    9 , 7.59 , 0.244 , 6.85 , 0.204
    10 , 7.60 , 0.254 , 6.74 , 0.230
    12 , 7.60 , 0.311 , 6.87 , 0.216
    14 , 7.71 , 0.250 , 7.30 , 0.192
    15 , 7.59 , 0.235 , 6.77 , 0.433
    18 , 7.57 , 0.291 , 6.90 , 0.227
    24 , 7.73 , 0.335 , 6.59 , 0.208
    25 , 7.46 , 0.275 , 7.14 , 0.222
    26 , 7.64 , 0.243 , 7.05 , 0.146
    32 , 7.70 , 0.303 , 6.94 , 0.282
    33 , 7.74 , 0.258 , 6.89 , 0.265
    35 , 7.62 , 0.258 , 7.10 , 0.222
    37 , 7.72 , 0.250 , 6.73 , 0.241
    38 , 7.62 , 0.230 , 6.93 , 0.215
    43 , 7.41 , 0.253 , 7.20 , 0.246
    44 , 7.26 , 0.230 , 7.05 , 0.251
    46 , 7.70 , 0.322 , 7.13 , 0.212
    47 , 7.47 , 0.316 , 7.40 , 0.263
    49 , 7.74 , 0.379 , 7.32 , 0.271
    53 , 7.63 , 0.211 , 6.48 , 0.171
    54 , 7.58 , 0.239 , 7.22 , 0.262
    55 , 7.41 , 0.249 , 6.51 , 0.221
    58 , 7.62 , 0.276 , 6.90 , 0.231
    60 , 7.64 , 0.226 , 6.62 , 0.196
    86 , 7.45 , 0.232 , 6.69 , 0.188
    87 , 7.42 , 0.246 , 7.11 , 0.227

