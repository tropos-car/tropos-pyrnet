## Calibration for IOPRAO 2014
```
processed with pyrnet-0.2.16
```

The PyrNet was setup for calibration at the DWD Lindenberg facility and the Falkenberg field from 2014-06-02 to 2014-07-18. Cross-calibration is done versus reference observations from the TROPOS MObile RaDiation ObseRvatory (MORDOR) station and BSRN measurement station at Lindeberg ([Wacker & Behrens 2022](https://doi.pangaea.de/10.1594/PANGAEA.946382)).

As PyrNet stations are not clusterd on a sigle facility, but several kilometers apart, only clear sky screened reference data is used for calibration, as the sun position differences are negligible.

### Imports


```python
#|dropcode
import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import jstyleson as json

from pvlib import clearsky
from pvlib.location import Location

from pyrnet import pyrnet
```

### Prepare PyrNet data
For calibration preparation the PyrNet data is processed to level l1b using a calibration factor of **7 (uV W-1 m2)** for all pyranometers with the ```pyrnet process l1b``` tool. This is done to unify the conversion to sensor voltage during calibration and not run into valid_range limits for netcdf encoding. Here we generate the *calibration.json* file for the processing to l1b:   


```python
box_numbers = np.arange(1,101)
calibrations = {f"{bn:03d}":[7,7] for bn in box_numbers}
calibjson = {"2000-01-01": calibrations}
with open("pyrnet_calib_prep.json","w") as txt:
    json.dump(calibjson, txt)
```

Within *pyrnet_config.json*:
```
{"file_calibration" : "pyrnet_calib_prep.json"}
```

**Workflow for preparation**
1. Prepare *pyrnet_config_calibration_prep.json* with contributors metadata and the dummy calibration config file.
1. ```$ pyrnet process l1a -c pyrnet_config.json raw_data/*.bin l1a/```
1. ```$ pyrnet process l1b_network -c pyrnet_config.json l1a/*.nc l1b_network/```

### Configuration
Set local data paths and lookup metadata.


```python
pf_mordor = "mordor/{date:%Y/%m/%Y-%m-%d}_Radiation.dat"
pf_bsrn = "bsrn/{date:%Y-%m}_bsrn.tab"
pf_pyrnet = "l1b_network/{date:%Y-%m-%d}_P1D_pyrnet_ioprao_n000l1bf1s.c01.nc"
dates = pd.date_range("2014-06-04","2014-07-18")
# period with lots of clear sky situations
dates = pd.date_range("2014-06-06","2014-06-09")
stations = np.arange(1,101)

loc = Location(52.21, 14.122, altitude=125) # Lindenberg

# lookup which box contains actually a pyranometer/ extra pyranometer
mainmask = [] 
for box in stations:
    _, serials, _, _ = pyrnet.meta_lookup(dates[0],box=box)
    mainmask.append( True if len(serials[0])>0 else False )
```

#### Load reference MORDOR data


```python
#|dropcode
#|dropout
new = True
for i,date in enumerate(dates):
    fname = pf_mordor.format(date=date)
    if not os.path.exists(fname):
        continue
    df = pd.read_csv(
        fname,
        header=0,
        skiprows=[0,2,3],
        date_format="ISO8601",
        na_values=["NAN"],
        parse_dates=[0],
        index_col=0
    )
    dst = df.to_xarray().rename({"TIMESTAMP":"time"})

    # drop not needed variables
    keep_vars = ['TP2_Wm2'] # global shortwave irradiance
    drop_vars = [v for v in dst if v not in keep_vars]
    dst = dst.drop_vars(drop_vars)
    dst = dst.resample(time="1min").mean(skipna=True)

    cs = loc.get_clearsky(pd.to_datetime(dst.time.values),model='simplified_solis')
    try:
        cs_mask = clearsky.detect_clearsky(
            dst['TP2_Wm2'].values,
            cs['ghi'],
            times=pd.to_datetime(dst.time.values)
        )    
    except:
        cs_mask = np.zeros(dst.time.size).astype(bool)

    dst = dst.assign({"cs_mask":("time", cs_mask)})
    
    # merge
    if new:
        ds = dst.copy()
        new = False
    else:
        ds = xr.concat((ds,dst),dim='time', data_vars='minimal', coords='minimal', compat='override')

mordor = ds.copy()
mordor = mordor.drop_duplicates("time", keep="last")
mordor
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body[data-theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block !important;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-index-preview {
  grid-column: 2 / 5;
  color: var(--xr-font-color2);
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data,
.xr-index-data-in:checked ~ .xr-index-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-index-name div,
.xr-index-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2,
.xr-no-icon {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:  (time: 5719)
Coordinates:
  * time     (time) datetime64[ns] 2014-06-06 ... 2014-06-09T23:18:00
Data variables:
    TP2_Wm2  (time) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... nan nan nan nan nan 0.0
    cs_mask  (time) bool False False False False ... False False False False</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-5ff5dcce-c90b-4684-9fdb-6b471b153a77' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-5ff5dcce-c90b-4684-9fdb-6b471b153a77' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 5719</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-58a54dd6-9dc9-454d-b321-873466fc78e5' class='xr-section-summary-in' type='checkbox'  checked><label for='section-58a54dd6-9dc9-454d-b321-873466fc78e5' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2014-06-06 ... 2014-06-09T23:18:00</div><input id='attrs-e037268d-7789-455a-8c8d-f709b528e0c4' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-e037268d-7789-455a-8c8d-f709b528e0c4' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-66f94d2f-d1fe-44df-8db1-c29afd3fb827' class='xr-var-data-in' type='checkbox'><label for='data-66f94d2f-d1fe-44df-8db1-c29afd3fb827' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2014-06-06T00:00:00.000000000&#x27;, &#x27;2014-06-06T00:01:00.000000000&#x27;,
       &#x27;2014-06-06T00:02:00.000000000&#x27;, ..., &#x27;2014-06-09T23:16:00.000000000&#x27;,
       &#x27;2014-06-09T23:17:00.000000000&#x27;, &#x27;2014-06-09T23:18:00.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-1bfb458d-d9ce-40e9-9fa7-d95ce2318876' class='xr-section-summary-in' type='checkbox'  checked><label for='section-1bfb458d-d9ce-40e9-9fa7-d95ce2318876' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>TP2_Wm2</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0 0.0 0.0 0.0 ... nan nan nan 0.0</div><input id='attrs-f35458fc-7ef8-4d35-abcd-42a5bd20c27e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-f35458fc-7ef8-4d35-abcd-42a5bd20c27e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-ff57c7a9-ed44-4235-89ed-94fb24a3776e' class='xr-var-data-in' type='checkbox'><label for='data-ff57c7a9-ed44-4235-89ed-94fb24a3776e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 0.,  0.,  0., ..., nan, nan,  0.])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>cs_mask</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>bool</div><div class='xr-var-preview xr-preview'>False False False ... False False</div><input id='attrs-910743db-4c70-4d6a-a48d-0dbe84937bc2' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-910743db-4c70-4d6a-a48d-0dbe84937bc2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-4ec9f8da-0000-4893-b621-5e791d13ed9c' class='xr-var-data-in' type='checkbox'><label for='data-4ec9f8da-0000-4893-b621-5e791d13ed9c' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([False, False, False, ..., False, False, False])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-db00bc68-40d7-4b7a-8905-81065a81bc21' class='xr-section-summary-in' type='checkbox'  ><label for='section-db00bc68-40d7-4b7a-8905-81065a81bc21' class='xr-section-summary' >Indexes: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-9e2ef95b-50e4-41d5-bb02-4e7a789a82a6' class='xr-index-data-in' type='checkbox'/><label for='index-9e2ef95b-50e4-41d5-bb02-4e7a789a82a6' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2014-06-06 00:00:00&#x27;, &#x27;2014-06-06 00:01:00&#x27;,
               &#x27;2014-06-06 00:02:00&#x27;, &#x27;2014-06-06 00:03:00&#x27;,
               &#x27;2014-06-06 00:04:00&#x27;, &#x27;2014-06-06 00:05:00&#x27;,
               &#x27;2014-06-06 00:06:00&#x27;, &#x27;2014-06-06 00:07:00&#x27;,
               &#x27;2014-06-06 00:08:00&#x27;, &#x27;2014-06-06 00:09:00&#x27;,
               ...
               &#x27;2014-06-09 23:09:00&#x27;, &#x27;2014-06-09 23:10:00&#x27;,
               &#x27;2014-06-09 23:11:00&#x27;, &#x27;2014-06-09 23:12:00&#x27;,
               &#x27;2014-06-09 23:13:00&#x27;, &#x27;2014-06-09 23:14:00&#x27;,
               &#x27;2014-06-09 23:15:00&#x27;, &#x27;2014-06-09 23:16:00&#x27;,
               &#x27;2014-06-09 23:17:00&#x27;, &#x27;2014-06-09 23:18:00&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=5719, freq=&#x27;T&#x27;))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-0973e383-bfeb-4729-b876-438d8796be5a' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-0973e383-bfeb-4729-b876-438d8796be5a' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
fig,ax = plt.subplots(1,1)
ax.plot(mordor.time,mordor.TP2_Wm2)
ax.plot(mordor.time[mordor.cs_mask],mordor.TP2_Wm2[mordor.cs_mask],ls='',marker='.')
```




    [<matplotlib.lines.Line2D at 0x7efd5e626050>]




    
![png](calibration_ioprao_output_10_1.png)
    


## Load BSRN data


```python
umonth = np.unique(dates.values.astype("datetime64[M]"))

new = True
for month in umonth:
    fname = pf_bsrn.format(date=pd.to_datetime(month))
    if not os.path.exists(fname):
        continue
    df = pd.read_csv(
        fname,
        sep='\s+',
        header=None,
        skiprows=34,
        date_format="ISO8601",
        parse_dates=[0],
        index_col=0,
        names=["time","SWD"],
        usecols=[0,2]
    )
    dst = df.to_xarray()

    dst.SWD.values[dst.SWD.values<0] = 0
    dst.SWD.values = dst.SWD.values.astype(float)

    cs = loc.get_clearsky(pd.to_datetime(dst.time.values),model='simplified_solis')
    
    try:
        cs_mask = clearsky.detect_clearsky(
            dst['SWD'].values,
            cs['ghi'].values,
            times=pd.to_datetime(dst.time.values),
            window_length=60
        ) 
    except:
        cs_mask = np.zeros(dst.time.size).astype(bool)

    dst = dst.assign({"cs_mask":("time", cs_mask)})
    
    # merge
    if new:
        ds = dst.copy()
        new = False
    else:
        ds = xr.concat((ds,dst),dim='time', data_vars='minimal', coords='minimal', compat='override')


bsrn = ds.copy()
bsrn = bsrn.drop_duplicates("time", keep="last")
bsrn
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body[data-theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block !important;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-index-preview {
  grid-column: 2 / 5;
  color: var(--xr-font-color2);
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data,
.xr-index-data-in:checked ~ .xr-index-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-index-name div,
.xr-index-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2,
.xr-no-icon {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:  (time: 43200)
Coordinates:
  * time     (time) datetime64[ns] 2014-06-01 ... 2014-06-30T23:59:00
Data variables:
    SWD      (time) float64 1.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0 0.0
    cs_mask  (time) bool False False False False ... False False False False</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-f87bee71-39d7-45a3-a489-af2854f5e0f9' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-f87bee71-39d7-45a3-a489-af2854f5e0f9' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 43200</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-2de6acb7-6dcc-42cc-8bbf-5fb7c877659a' class='xr-section-summary-in' type='checkbox'  checked><label for='section-2de6acb7-6dcc-42cc-8bbf-5fb7c877659a' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2014-06-01 ... 2014-06-30T23:59:00</div><input id='attrs-72bab6f0-07d3-4f86-8a98-befce89c7e50' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-72bab6f0-07d3-4f86-8a98-befce89c7e50' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-63911c65-6667-4936-80e1-6d1ef597a904' class='xr-var-data-in' type='checkbox'><label for='data-63911c65-6667-4936-80e1-6d1ef597a904' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2014-06-01T00:00:00.000000000&#x27;, &#x27;2014-06-01T00:01:00.000000000&#x27;,
       &#x27;2014-06-01T00:02:00.000000000&#x27;, ..., &#x27;2014-06-30T23:57:00.000000000&#x27;,
       &#x27;2014-06-30T23:58:00.000000000&#x27;, &#x27;2014-06-30T23:59:00.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-f9faef45-332c-435f-86f1-233a190e4507' class='xr-section-summary-in' type='checkbox'  checked><label for='section-f9faef45-332c-435f-86f1-233a190e4507' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>SWD</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0</div><input id='attrs-f5a6a482-4a86-429a-b448-17164000b3f6' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-f5a6a482-4a86-429a-b448-17164000b3f6' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-e3ecdf9d-1b93-4a84-9714-6119ebbe8d85' class='xr-var-data-in' type='checkbox'><label for='data-e3ecdf9d-1b93-4a84-9714-6119ebbe8d85' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([1., 0., 0., ..., 0., 0., 0.])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>cs_mask</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>bool</div><div class='xr-var-preview xr-preview'>False False False ... False False</div><input id='attrs-c07aeea0-e0da-4382-b7aa-d15f4f56d5f8' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-c07aeea0-e0da-4382-b7aa-d15f4f56d5f8' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-230a355b-1639-4f18-9a8a-cfb51c5ced63' class='xr-var-data-in' type='checkbox'><label for='data-230a355b-1639-4f18-9a8a-cfb51c5ced63' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([False, False, False, ..., False, False, False])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-edc52274-900b-42cc-8bd0-cd5449c6003b' class='xr-section-summary-in' type='checkbox'  ><label for='section-edc52274-900b-42cc-8bd0-cd5449c6003b' class='xr-section-summary' >Indexes: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-e7b03aaf-f49a-4d75-b29c-196d626c5fbb' class='xr-index-data-in' type='checkbox'/><label for='index-e7b03aaf-f49a-4d75-b29c-196d626c5fbb' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2014-06-01 00:00:00&#x27;, &#x27;2014-06-01 00:01:00&#x27;,
               &#x27;2014-06-01 00:02:00&#x27;, &#x27;2014-06-01 00:03:00&#x27;,
               &#x27;2014-06-01 00:04:00&#x27;, &#x27;2014-06-01 00:05:00&#x27;,
               &#x27;2014-06-01 00:06:00&#x27;, &#x27;2014-06-01 00:07:00&#x27;,
               &#x27;2014-06-01 00:08:00&#x27;, &#x27;2014-06-01 00:09:00&#x27;,
               ...
               &#x27;2014-06-30 23:50:00&#x27;, &#x27;2014-06-30 23:51:00&#x27;,
               &#x27;2014-06-30 23:52:00&#x27;, &#x27;2014-06-30 23:53:00&#x27;,
               &#x27;2014-06-30 23:54:00&#x27;, &#x27;2014-06-30 23:55:00&#x27;,
               &#x27;2014-06-30 23:56:00&#x27;, &#x27;2014-06-30 23:57:00&#x27;,
               &#x27;2014-06-30 23:58:00&#x27;, &#x27;2014-06-30 23:59:00&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=43200, freq=None))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-144465d6-6309-4f26-b172-5cb3b4b85348' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-144465d6-6309-4f26-b172-5cb3b4b85348' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
bsrnr = bsrn.interp_like(mordor)

fig,ax = plt.subplots(1,1)
ax.set_title("Compare BSRN and MORDOR")
ax.scatter(bsrnr.SWD.values, mordor.TP2_Wm2.values, alpha=0.1 ,color='k')
ax.axline((0,0),slope=1,c='r',ls='--')
ax.set_ylabel("MORDOR GHI (Wm-2)")
ax.set_xlabel("BSRN GHI (Wm-2)")

ax.grid(True)
```


    
![png](calibration_ioprao_output_13_0.png)
    


#### Load PyrNet Data


```python
#|dropcode
#|dropout
for i,date in enumerate(dates):
    # read from thredds server
    dst = xr.open_dataset(pf_pyrnet.format(date=date))
    
    # drop not needed variables
    keep_vars = ['ghi','szen']
    drop_vars = [v for v in dst if v not in keep_vars]
    dst = dst.drop_vars(drop_vars)

    # unify time and station dimension to speed up merging
    date = dst.time.values[0].astype("datetime64[D]")
    timeidx = pd.date_range(date, date + np.timedelta64(1, 'D'), freq='1s', inclusive='left')
    dst = dst.interp(time=timeidx)
    dst = dst.reindex({"station": stations})

    dst.ghi.values = dst.ghi.values * 7 * 1e-6
    dst = dst.where(dst.szen<80, drop=True)
    dst.ghi.values = dst.ghi.where(dst.ghi>0.033/300.).values

    dst = dst.resample(time="1min").mean(skipna=True)
    
    # merge
    if i == 0:
        ds = dst.copy()
    else:
        ds = xr.concat((ds,dst),dim='time', data_vars='minimal', coords='minimal', compat='override')
    
pyr = ds.copy()
pyr
```




<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body[data-theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block !important;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-index-preview {
  grid-column: 2 / 5;
  color: var(--xr-font-color2);
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data,
.xr-index-data-in:checked ~ .xr-index-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-index-name div,
.xr-index-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data,
.xr-index-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2,
.xr-no-icon {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:          (time: 3319, station: 43, maintenancetime: 1)
Coordinates:
  * station          (station) int64 1 5 6 8 15 16 20 ... 90 91 92 93 95 98 99
  * maintenancetime  (maintenancetime) datetime64[ns] 2014-06-11T23:59:59
  * time             (time) datetime64[ns] 2014-06-06T04:08:00 ... 2014-06-09...
Data variables:
    ghi              (time, station) float64 0.0003242 0.0003215 ... nan nan
    szen             (time, station) float64 79.99 79.99 79.99 ... 79.99 nan nan
Attributes: (12/31)
    Conventions:               CF-1.10, ACDD-1.3
    title:                     TROPOS pyranometer network (PyrNet) observatio...
    history:                   2024-11-12T10:03:39: Merged level l1b by pyrne...
    institution:               Leibniz Institute for Tropospheric Research (T...
    source:                    TROPOS pyranometer network (PyrNet)
    references:                https://doi.org/10.5194/amt-9-1153-2016
    ...                        ...
    geospatial_lon_units:      degE
    time_coverage_start:       2014-06-06T00:00:00
    time_coverage_end:         2014-06-06T23:59:59
    time_coverage_duration:    P0DT23H59M59S
    time_coverage_resolution:  P0DT0H0M1S
    site:                      [&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;...</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-42694a30-5c1d-44f4-9b7e-da4d410c604c' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-42694a30-5c1d-44f4-9b7e-da4d410c604c' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 3319</li><li><span class='xr-has-index'>station</span>: 43</li><li><span class='xr-has-index'>maintenancetime</span>: 1</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-a591f432-b3b8-4562-9f0b-4571129eeddd' class='xr-section-summary-in' type='checkbox'  checked><label for='section-a591f432-b3b8-4562-9f0b-4571129eeddd' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>station</span></div><div class='xr-var-dims'>(station)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 5 6 8 15 16 ... 91 92 93 95 98 99</div><input id='attrs-9e0485b0-9fdd-4280-9208-d8d515cc35d5' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-9e0485b0-9fdd-4280-9208-d8d515cc35d5' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-0fa3da03-eb21-4a9d-9f89-3a19d56f707b' class='xr-var-data-in' type='checkbox'><label for='data-0fa3da03-eb21-4a9d-9f89-3a19d56f707b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 1,  5,  6,  8, 15, 16, 20, 21, 22, 30, 37, 41, 42, 43, 44, 45, 46, 47,
       49, 50, 53, 54, 55, 56, 57, 61, 64, 68, 71, 74, 75, 80, 82, 86, 87, 88,
       90, 91, 92, 93, 95, 98, 99])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>maintenancetime</span></div><div class='xr-var-dims'>(maintenancetime)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2014-06-11T23:59:59</div><input id='attrs-5da39f38-47a5-4f9f-9ea0-9aca645bba5c' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-5da39f38-47a5-4f9f-9ea0-9aca645bba5c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-a1ffd7e8-037c-4c47-b026-156ec05247ae' class='xr-var-data-in' type='checkbox'><label for='data-a1ffd7e8-037c-4c47-b026-156ec05247ae' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2014-06-11T23:59:59.000000000&#x27;], dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2014-06-06T04:08:00 ... 2014-06-...</div><input id='attrs-2b9c3307-bad8-49a0-b3bf-75550a7ab62c' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-2b9c3307-bad8-49a0-b3bf-75550a7ab62c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-bfb0a3ff-331e-474d-b92b-66da495eeabb' class='xr-var-data-in' type='checkbox'><label for='data-bfb0a3ff-331e-474d-b92b-66da495eeabb' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2014-06-06T04:08:00.000000000&#x27;, &#x27;2014-06-06T04:09:00.000000000&#x27;,
       &#x27;2014-06-06T04:10:00.000000000&#x27;, ..., &#x27;2014-06-09T17:56:00.000000000&#x27;,
       &#x27;2014-06-09T17:57:00.000000000&#x27;, &#x27;2014-06-09T17:58:00.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-4848ee5d-2c5b-4f27-95bd-97b2804679ae' class='xr-section-summary-in' type='checkbox'  checked><label for='section-4848ee5d-2c5b-4f27-95bd-97b2804679ae' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>ghi</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0003242 0.0003215 ... nan nan</div><input id='attrs-1fa49c21-b2ad-4679-9ffe-109a91bd9477' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-1fa49c21-b2ad-4679-9ffe-109a91bd9477' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-90a2bb14-50a5-4039-9bc7-531c0bd3dd95' class='xr-var-data-in' type='checkbox'><label for='data-90a2bb14-50a5-4039-9bc7-531c0bd3dd95' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>W m-2</dd><dt><span>long_name :</span></dt><dd>downwelling shortwave flux</dd><dt><span>standard_name :</span></dt><dd>downwelling_shortwave_flux_in_air</dd><dt><span>valid_range :</span></dt><dd>[    0 60000]</dd><dt><span>ancillary_variables :</span></dt><dd>maintenance_flag_ghi qc_flag_ghi</dd><dt><span>serial :</span></dt><dd>[&#x27;S12128.001&#x27;, &#x27;S12128.005&#x27;, &#x27;S12078.061&#x27;, &#x27;S12128.008&#x27;, &#x27;S12128.015&#x27;, &#x27;S12128.027&#x27;, &#x27;S12128.020&#x27;, &#x27;S12128.021&#x27;, &#x27;S12128.022&#x27;, &#x27;S12128.030&#x27;, &#x27;S12128.037&#x27;, &#x27;S12128.036&#x27;, &#x27;S12128.042&#x27;, &#x27;S12128.043&#x27;, &#x27;S12128.044&#x27;, &#x27;S12128.045&#x27;, &#x27;S12128.046&#x27;, &#x27;S12128.047&#x27;, &#x27;S12128.049&#x27;, &#x27;S12128.050&#x27;, &#x27;S12137.003&#x27;, &#x27;S12137.004&#x27;, &#x27;S12137.005&#x27;, &#x27;S12137.006&#x27;, &#x27;S12137.007&#x27;, &#x27;S12137.011&#x27;, &#x27;S12137.014&#x27;, &#x27;S12137.018&#x27;, &#x27;S12137.021&#x27;, &#x27;S12137.024&#x27;, &#x27;S12137.025&#x27;, &#x27;S12137.030&#x27;, &#x27;S12137.032&#x27;, &#x27;S12137.036&#x27;, &#x27;S12137.037&#x27;, &#x27;S12137.038&#x27;, &#x27;S12137.040&#x27;, &#x27;S12137.041&#x27;, &#x27;S12137.042&#x27;, &#x27;S12137.043&#x27;, &#x27;S12137.045&#x27;, &#x27;S12137.048&#x27;, &#x27;S12137.049&#x27;]</dd><dt><span>calibration_Cabsolute :</span></dt><dd>[142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714]</dd><dt><span>calibration_Ccoscorr :</span></dt><dd>1.0</dd><dt><span>calibration_function :</span></dt><dd>flux (W m-2) = flux (V) * Cabsolute (W m-2 V-1) * Ccoscorr(mua)</dd><dt><span>vangle :</span></dt><dd>[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0]</dd><dt><span>hangle :</span></dt><dd>[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0]</dd></dl></div><div class='xr-var-data'><pre>array([[0.00032423, 0.00032148, 0.00030234, ..., 0.00034445, 0.0003227 ,
        0.00034237],
       [0.00033008, 0.00032603, 0.00030423, ..., 0.00035312, 0.000325  ,
        0.000347  ],
       [0.00032804, 0.00032221, 0.0003091 , ..., 0.00036259, 0.00033331,
        0.00035446],
       ...,
       [0.00048062, 0.00047669, 0.00040849, ..., 0.00043054, 0.00044911,
        0.00045594],
       [0.00049336, 0.00048807, 0.00042757, ..., 0.00043078, 0.00046371,
        0.00047683],
       [       nan,        nan,        nan, ..., 0.00043715,        nan,
               nan]])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>szen</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>79.99 79.99 79.99 ... 79.99 nan nan</div><input id='attrs-2fef65fa-d0e4-4222-a39e-fa3ce9199697' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-2fef65fa-d0e4-4222-a39e-fa3ce9199697' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-de9cfe02-f77f-4323-a9ad-3c167e78b4f4' class='xr-var-data-in' type='checkbox'><label for='data-de9cfe02-f77f-4323-a9ad-3c167e78b4f4' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>standard_name :</span></dt><dd>solar_zenith_angle</dd><dt><span>units :</span></dt><dd>degree</dd><dt><span>valid_range :</span></dt><dd>[    0 36000]</dd></dl></div><div class='xr-var-data'><pre>array([[79.99249649, 79.99249649, 79.99166361, ..., 79.98192127,
        79.99166361, 79.99099731],
       [79.91716487, 79.91749827, 79.91449814, ..., 79.89683151,
        79.91474813, 79.9154981 ],
       [79.77708155, 79.77741496, 79.77449824, ..., 79.75683136,
        79.77449824, 79.77524821],
       ...,
       [79.79641482, 79.79616483, 79.80174815, ..., 79.78299828,
        79.80174815, 79.79866486],
       [79.93219101, 79.93201554, 79.93407228, ..., 79.92224808,
        79.93407228, 79.93348026],
       [        nan,         nan,         nan, ..., 79.99499512,
                nan,         nan]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-1b198b4b-03f7-4792-8b9f-58f92d439f05' class='xr-section-summary-in' type='checkbox'  ><label for='section-1b198b4b-03f7-4792-8b9f-58f92d439f05' class='xr-section-summary' >Indexes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>station</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-ff6ce919-6819-435e-a6ae-ab52e6efecfe' class='xr-index-data-in' type='checkbox'/><label for='index-ff6ce919-6819-435e-a6ae-ab52e6efecfe' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([ 1,  5,  6,  8, 15, 16, 20, 21, 22, 30, 37, 41, 42, 43, 44, 45, 46, 47,
       49, 50, 53, 54, 55, 56, 57, 61, 64, 68, 71, 74, 75, 80, 82, 86, 87, 88,
       90, 91, 92, 93, 95, 98, 99],
      dtype=&#x27;int64&#x27;, name=&#x27;station&#x27;))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>maintenancetime</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-892ffaec-a52f-41a9-9737-09d8d35eec80' class='xr-index-data-in' type='checkbox'/><label for='index-892ffaec-a52f-41a9-9737-09d8d35eec80' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2014-06-11 23:59:59&#x27;], dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;maintenancetime&#x27;, freq=None))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-ea2da4d2-47bd-4fa1-9a2f-430e8ccb28b1' class='xr-index-data-in' type='checkbox'/><label for='index-ea2da4d2-47bd-4fa1-9a2f-430e8ccb28b1' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2014-06-06 04:08:00&#x27;, &#x27;2014-06-06 04:09:00&#x27;,
               &#x27;2014-06-06 04:10:00&#x27;, &#x27;2014-06-06 04:11:00&#x27;,
               &#x27;2014-06-06 04:12:00&#x27;, &#x27;2014-06-06 04:13:00&#x27;,
               &#x27;2014-06-06 04:14:00&#x27;, &#x27;2014-06-06 04:15:00&#x27;,
               &#x27;2014-06-06 04:16:00&#x27;, &#x27;2014-06-06 04:17:00&#x27;,
               ...
               &#x27;2014-06-09 17:49:00&#x27;, &#x27;2014-06-09 17:50:00&#x27;,
               &#x27;2014-06-09 17:51:00&#x27;, &#x27;2014-06-09 17:52:00&#x27;,
               &#x27;2014-06-09 17:53:00&#x27;, &#x27;2014-06-09 17:54:00&#x27;,
               &#x27;2014-06-09 17:55:00&#x27;, &#x27;2014-06-09 17:56:00&#x27;,
               &#x27;2014-06-09 17:57:00&#x27;, &#x27;2014-06-09 17:58:00&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=3319, freq=None))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-3093d086-707b-4eaa-a23c-151463e42c20' class='xr-section-summary-in' type='checkbox'  ><label for='section-3093d086-707b-4eaa-a23c-151463e42c20' class='xr-section-summary' >Attributes: <span>(31)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>Conventions :</span></dt><dd>CF-1.10, ACDD-1.3</dd><dt><span>title :</span></dt><dd>TROPOS pyranometer network (PyrNet) observational data set</dd><dt><span>history :</span></dt><dd>2024-11-12T10:03:39: Merged level l1b by pyrnet version 0.2.16+9.gbd94801; </dd><dt><span>institution :</span></dt><dd>Leibniz Institute for Tropospheric Research (TROPOS)</dd><dt><span>source :</span></dt><dd>TROPOS pyranometer network (PyrNet)</dd><dt><span>references :</span></dt><dd>https://doi.org/10.5194/amt-9-1153-2016</dd><dt><span>Department :</span></dt><dd>Remote Sensing of Atmospheric Processes</dd><dt><span>Department_team :</span></dt><dd>Clouds, Aerosol and Radiation</dd><dt><span>Address :</span></dt><dd>Permoser Str. 15, 04318 Leipzig, Germany</dd><dt><span>Contact_person :</span></dt><dd>Andreas Macke and the clouds, aerosol and radiation team of the remote sensing department, mailto:andreas.macke@tropos.de</dd><dt><span>Contributor_name :</span></dt><dd>[&#x27;Leibniz-Institute for Tropospheric Research (TROPOS)&#x27;, &#x27;German Weather Service, Observatory Lindenberg (DWD)&#x27;, &#x27;Andreas Macke (TROPOS)&#x27;, &#x27;Hartwig Deneke (TROPOS)&#x27;, &#x27;Bomidi Lakshmi Madhavan (TROPOS)&#x27;, &#x27;Jonas Witthuhn (TROPOS)&#x27;]</dd><dt><span>Contributor_role :</span></dt><dd>[&#x27;supervision, logistics, setup, maintenance, tear-down, data-processing, site-planning&#x27;, &#x27;supervision, site-manager&#x27;, &#x27;supervision&#x27;, &#x27;supervision&#x27;, &#x27;site-planning, logistics, setup, tear-down, data-processing, quality-control&#x27;, &#x27;logistics, setup, tear-down, maintenance, data-processing, quality-control&#x27;]</dd><dt><span>Authors_software :</span></dt><dd>Hartwig Deneke, Jonas Witthuhn, mailto:deneke@tropos.de</dd><dt><span>Creator_name :</span></dt><dd>Jonas Witthuhn</dd><dt><span>Project :</span></dt><dd>HD(CP)2 short term observations, Intensive observation periode RAO (IOPRAO)</dd><dt><span>Standard_name_vocabulary :</span></dt><dd>CF Standard Name Table v81</dd><dt><span>License :</span></dt><dd>CC-BY-SA 3.0</dd><dt><span>processing_level :</span></dt><dd>l1b</dd><dt><span>product_version :</span></dt><dd>0.2.16+9.gbd94801</dd><dt><span>date_created :</span></dt><dd>2024-11-12T10:03:40</dd><dt><span>geospatial_lat_min :</span></dt><dd>52.1643722272286</dd><dt><span>geospatial_lat_max :</span></dt><dd>52.209451299476896</dd><dt><span>geospatial_lat_units :</span></dt><dd>degN</dd><dt><span>geospatial_lon_min :</span></dt><dd>14.114361177601001</dd><dt><span>geospatial_lon_max :</span></dt><dd>14.128804051166476</dd><dt><span>geospatial_lon_units :</span></dt><dd>degE</dd><dt><span>time_coverage_start :</span></dt><dd>2014-06-06T00:00:00</dd><dt><span>time_coverage_end :</span></dt><dd>2014-06-06T23:59:59</dd><dt><span>time_coverage_duration :</span></dt><dd>P0DT23H59M59S</dd><dt><span>time_coverage_resolution :</span></dt><dd>P0DT0H0M1S</dd><dt><span>site :</span></dt><dd>[&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;]</dd></dl></div></li></ul></div></div>



### Calibration
The calibration follows the [ISO 9847:1992 - Solar energy — Calibration of field pyranometers by comparison to a reference pyranometer](https://archive.org/details/gov.in.is.iso.9847.1992).
> TODO: Revise versus 2023 EU version.

Cloudy sky treatment is applied.

#### Step 1
Drop Night measures and low signal measures from pyranometer data. Since calibration without incoming radiation doesnt work.

This data is kept for calibration:

 * solar zenith angle < 80° ( as recommended in ISO 9847)
 * Measured Voltage > 0.033 V, e.g. ADC count is 0 or 1 of 1023 (drop the lowest ~1%)
 
Voltage measured ($V_m$) at the logger is the amplified Senor voltage ($V_S$) by a gain of 300.

$ V_m = 300 V_S$

As the uncalibrated flux measurements ($F_U$) are calibrated with a fixed factor of 7 uV W-1 m2:

$ V_s = 7*1e-6* F_U $


```python
# # Set flux values to nan if no pyranometer is installed.
# pyr.ghi.values = pyr.ghi.where(mainmask).values

```

#### Step 2
Interpolate reference to PyrNet samples and combine to a single Dataset


```python
# interpolate reference to PyrNet
mordor = mordor.interp(time=pyr.time).interpolate_na()
bsrn = bsrn.interp(time=pyr.time).interpolate_na()
```


```python
fig,ax = plt.subplots(1,1)
ax.plot(mordor.time,mordor.TP2_Wm2)
ax.plot(mordor.time[mordor.cs_mask],mordor.TP2_Wm2[mordor.cs_mask],ls='',marker='.')
ax.plot(bsrn.time,bsrn.SWD+100)
ax.plot(bsrn.time[bsrn.cs_mask],bsrn.SWD[bsrn.cs_mask]+100,ls='',marker='.')
```




    [<matplotlib.lines.Line2D at 0x7efd5e63b190>]




    
![png](calibration_ioprao_output_21_1.png)
    



```python
# Calibration datasets for main and extra pyranometer
Cds_main = xr.Dataset(
    data_vars={
        'reference2_Wm2': ('time', mordor['TP2_Wm2'].data[bsrn.cs_mask]),
        'reference_Wm2': ('time', bsrn['SWD'].data[bsrn.cs_mask]),
        'pyrnet_V': (('time','station'), pyr['ghi'].data[bsrn.cs_mask])
    },
    coords= {
        "time": pyr.time[bsrn.cs_mask],
        "station": pyr.station
    }
)
```

#### Step 3
Remove outliers from series using xarray grouping and apply function. The following functions removes outliers (deviation more than 2% according to ISO 9847) from a selected group. This step involves calculating calibration series and the integration of one hour intervals to smooth out high variable situation, which would break the calibration even when time synchronization is slightly off. Also this gets rid of some random shading events like birds / chimney / rods in line of sigth, which would affect calibration otherwise. We following ISO 9847 5.4.1.1 equation (2) here.


```python
def remove_outliers(x):
    """
    x is an xarray dataset containing these variables:
    coords: 'time' - datetime64
    'pyrnet_V' - array - voltage measures of pyranometer
    'reference_Wm2' - array - measured irradiance of reference
    """

    # calculate calibration series for single samples
    C = x['pyrnet_V'] / x['reference_Wm2']
    # integrated series 
    ix = x.integrate('time')
    M = ix['pyrnet_V'] / ix['reference_Wm2']
    
    while np.any(np.abs(C-M) > 0.02*M):
        #calculate as long there are outliers deviating more than 2 percent
        x = x.where(np.abs(C-M) < 0.02*M)
        C = x['pyrnet_V'] / x['reference_Wm2']
        #integrated series 
        ix = x.integrate('time')
        M = ix['pyrnet_V'] / ix['reference_Wm2']
        
    #return the reduced dataset x
    return x

# remove outliers
Cds_main = Cds_main.groupby('time.hour').apply(remove_outliers)

# hourly mean
Cds_main = Cds_main.resample(time="1h").mean(skipna=True)

```

#### Step 4
The series of measured voltage and irradiance is now without outliers. So we use equation 1 again to calculate from this reduced series the calibration factor for the instant samples.


```python
C_main = 1e6*Cds_main['pyrnet_V'] / Cds_main['reference_Wm2']
C_main2 = 1e6*Cds_main['pyrnet_V'] / Cds_main['reference2_Wm2']
```


```python
C_main.values[C_main.values<6]=np.nan
C_main.values[C_main.values>8]=np.nan
C_main2.values[C_main2.values<6]=np.nan
C_main2.values[C_main2.values>8]=np.nan
```

#### Step 5
We just found the Calibration factor to be the mean of the reduced calibration factor series and the uncertainty to be the standard deviation of this reduced series. Steo 3, 4 and 5 are done for every pyranometer seperate.


```python
C_main_mean = C_main.mean(dim='time',skipna=True)
C_main_std = C_main.std(dim='time',skipna=True)

C_main2_mean = C_main2.mean(dim='time',skipna=True)
C_main2_std = C_main2.std(dim='time',skipna=True)
```

### Results


```python
#|dropcode
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_title("Main Pyranometer")
ax.plot(C_main.time, C_main, ls ="", marker='.',c='k')
ax.plot(C_main2.time, C_main2, ls ="", marker='x',c='b')
ax.set_xlabel("Date")
ax.set_ylabel("Calibration factor (uV / Wm-2)")
ax.grid(True)
fig.show()

plt.figure()
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_title("Main Pyranometer")
ax.plot(pyr.szen.interp_like(C_main), C_main, ls ="", marker='.',c='k')
ax.plot(pyr.szen.interp_like(C_main2), C_main2, ls ="", marker='x',c='b')
ax.set_xlabel("solar zenith angle (deg)")
ax.set_ylabel("Calibration factor (uV / Wm-2)")
ax.grid(True)
fig.show()



```


    
![png](calibration_ioprao_output_31_0.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](calibration_ioprao_output_31_2.png)
    



```python
calibration_new = {}
calibration2_new = {}
print(f"Box:    vs BSRN       ,      vs MORDOR ")
for box in C_main_mean.station:
    Cm = float(C_main_mean.sel(station=box).values)
    Um = float(C_main_std.sel(station=box).values)
    Cm2 = float(C_main2_mean.sel(station=box).values)
    Um2 = float(C_main2_std.sel(station=box).values)
    calibration_new.update({
        f"{box:03d}": [np.round(Cm,2), None]
    })
    calibration2_new.update({
        f"{box:03d}": [np.round(Cm2,2), None]
    })
    print(f"{box:3d}: {Cm:.2f} +- {Um:.3f} , {Cm2:.2f} +- {Um2:.3f}")

calibjson = {"2014-06-01": calibration_new}
with open("pyrnet_calib_new_bsrn.json","w") as txt:
    json.dump(calibjson, txt)

calibjson = {"2014-06-01": calibration2_new}
with open("pyrnet_calib_new_mordor.json","w") as txt:
    json.dump(calibjson, txt)
```

    Box:    vs BSRN       ,      vs MORDOR 
      1: 7.44 +- 0.198 , 7.45 +- 0.201
      5: 7.39 +- 0.300 , 7.35 +- 0.266
      6: 6.71 +- 0.253 , 6.73 +- 0.269
      8: 7.57 +- 0.346 , 7.51 +- 0.318
     15: 7.44 +- 0.217 , 7.47 +- 0.216
     16: 7.63 +- 0.190 , 7.61 +- 0.183
     20: 7.56 +- 0.196 , 7.55 +- 0.200
     21: 7.34 +- 0.224 , 7.39 +- 0.268
     22: 7.40 +- 0.219 , 7.42 +- 0.234
     30: 7.54 +- 0.230 , 7.47 +- 0.160
     37: 7.53 +- 0.174 , 7.54 +- 0.180
     41: 7.61 +- 0.206 , 7.57 +- 0.180
     42: 7.42 +- 0.280 , 7.45 +- 0.282
     43: 7.28 +- 0.161 , 7.33 +- 0.218
     44: 7.20 +- 0.399 , 7.18 +- 0.400
     45: 7.44 +- 0.243 , 7.45 +- 0.248
     46: 7.52 +- 0.261 , 7.53 +- 0.254
     47: 7.41 +- 0.254 , 7.37 +- 0.217
     49: 7.51 +- 0.225 , 7.46 +- 0.196
     50: 7.58 +- 0.199 , 7.59 +- 0.195
     53: 7.40 +- 0.301 , 7.42 +- 0.312
     54: 7.36 +- 0.305 , 7.39 +- 0.323
     55: 7.16 +- 0.237 , 7.18 +- 0.242
     56: 6.75 +- 0.188 , 6.77 +- 0.205
     57: 6.71 +- 0.198 , 6.73 +- 0.202
     61: 7.31 +- 0.205 , 7.33 +- 0.224
     64: 7.22 +- 0.178 , 7.24 +- 0.202
     68: 6.80 +- 0.142 , 6.82 +- 0.155
     71: 7.18 +- 0.228 , 7.21 +- 0.234
     74: 7.38 +- 0.080 , 7.40 +- 0.098
     75: 6.64 +- 0.158 , 6.63 +- 0.208
     80: 6.90 +- 0.176 , 6.92 +- 0.187
     82: 7.42 +- 0.187 , 7.43 +- 0.175
     86: 7.37 +- 0.157 , 7.40 +- 0.186
     87: 7.32 +- 0.205 , 7.34 +- 0.224
     88: 6.90 +- 0.205 , 6.92 +- 0.183
     90: 7.36 +- 0.218 , 7.38 +- 0.204
     91: 7.40 +- 0.244 , 7.45 +- 0.273
     92: 7.08 +- 0.280 , 7.16 +- 0.196
     93: 7.50 +- 0.161 , 7.52 +- 0.189
     95: 7.44 +- 0.253 , 7.46 +- 0.278
     98: 7.13 +- 0.240 , 7.15 +- 0.248
     99: 7.38 +- 0.193 , 7.40 +- 0.207

