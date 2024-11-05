## Calibration for MelCol 2015
```
processed with pyrnet-0.2.16
```

The PyrNet was setup for calibration in a dense array on the Melpitz measurement field from 2015-05-06 to 2015-05-11. Cross-calibration is done versus reference observations from the TROPOS MObile RaDiation ObseRvatory (MORDOR) station.

### Imports


```python
#|dropcode
from IPython.display import display, Latex
import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import jstyleson as json

from scipy.optimize import differential_evolution

from pvlib import clearsky
from pvlib.location import Location

import pyrnet.pyrnet
import pyrnet.utils
```

### Prepare PyrNet data
For calibration preparation the PyrNet data is processed to level l1b using a calibration factor of **7 (uV W-1 m2)** for all pyranometers with the ```pyrnet process l1b``` tool. This is done to unify the conversion to sensor voltage during calibration and not run into valid_range limits for netcdf encoding. 

> Before running this notebook, new absolute calibration factors have to be determined with calibration_melcol.ipynb

### Configuration
Set local data paths and lookup metadata.


```python
pf_mordor = "mordor/{date:%Y-%m-%d}_Radiation.dat"
pf_pyrnet = "l1b_network/{date:%Y-%m-%d}_P1D_pyrnet_melcol_n000l1bf1s.c01.nc"
dates = pd.date_range("2015-05-06","2015-05-11")

loc = Location(51.525175, 12.91648, altitude=90) # Melpitz

stations = np.arange(1,101)
# lookup which box contains actually a pyranometer/ extra pyranometer
mainmask = [] 
for box in stations:
    _, serials, _, _ = pyrnet.pyrnet.meta_lookup(dates[0],box=box)
    mainmask.append( True if len(serials[0])>0 else False )
```

#### Load reference MORDOR data
The reference data of MORDOR is loaded, and clearsky is detected on daily basis using solis_simple clearsky model and the pvlib.clearsky.detect_clearsky function.


```python
#|dropcode
#|dropout
for i,date in enumerate(dates):
    fname = pf_mordor.format(date=date)
    df = pd.read_csv(
        fname,
        header=0,
        skiprows=[0,2,3],
        date_format="ISO8601",
        parse_dates=[0],
        index_col=0
    )
    dst = df.to_xarray().rename({"TIMESTAMP":"time"})

    # drop not needed variables
    keep_vars = ['TP2_Wm2'] # GHI,DHI,DNI
    drop_vars = [v for v in dst if v not in keep_vars]
    dst = dst.drop_vars(drop_vars)
    dst = dst.drop_duplicates("time",keep="last")
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
    if i == 0:
        ds = dst.copy()
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
Dimensions:  (time: 8640)
Coordinates:
  * time     (time) datetime64[ns] 2015-05-06 ... 2015-05-11T23:59:00
Data variables:
    TP2_Wm2  (time) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0 0.0
    cs_mask  (time) bool False False False False ... False False False False</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-4849bc83-7dcf-4ee6-b7b5-3bd6607cee32' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-4849bc83-7dcf-4ee6-b7b5-3bd6607cee32' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 8640</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-ae27a2a0-3252-4400-96f5-e10edd9692d4' class='xr-section-summary-in' type='checkbox'  checked><label for='section-ae27a2a0-3252-4400-96f5-e10edd9692d4' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-06 ... 2015-05-11T23:59:00</div><input id='attrs-ecdbb1a6-3f45-44b3-af64-5bed3d9449c0' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-ecdbb1a6-3f45-44b3-af64-5bed3d9449c0' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-6dc7a4ab-5d1b-4ad8-92b5-56186d1ea416' class='xr-var-data-in' type='checkbox'><label for='data-6dc7a4ab-5d1b-4ad8-92b5-56186d1ea416' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-06T00:00:00.000000000&#x27;, &#x27;2015-05-06T00:01:00.000000000&#x27;,
       &#x27;2015-05-06T00:02:00.000000000&#x27;, ..., &#x27;2015-05-11T23:57:00.000000000&#x27;,
       &#x27;2015-05-11T23:58:00.000000000&#x27;, &#x27;2015-05-11T23:59:00.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-fb2ca093-bac7-49ba-bb5b-52ca0ed051fb' class='xr-section-summary-in' type='checkbox'  checked><label for='section-fb2ca093-bac7-49ba-bb5b-52ca0ed051fb' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>TP2_Wm2</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0</div><input id='attrs-d6038e6c-2b1b-47b2-b4f7-110601da9d7d' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-d6038e6c-2b1b-47b2-b4f7-110601da9d7d' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-26a40d33-f870-4b85-9dc7-cbc1e0467618' class='xr-var-data-in' type='checkbox'><label for='data-26a40d33-f870-4b85-9dc7-cbc1e0467618' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([0., 0., 0., ..., 0., 0., 0.])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>cs_mask</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>bool</div><div class='xr-var-preview xr-preview'>False False False ... False False</div><input id='attrs-e40d7e7a-cb4d-4adc-94eb-5f52f053060c' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-e40d7e7a-cb4d-4adc-94eb-5f52f053060c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-d8bd0f4b-e238-4005-8eb5-ce91d5b484bd' class='xr-var-data-in' type='checkbox'><label for='data-d8bd0f4b-e238-4005-8eb5-ce91d5b484bd' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([False, False, False, ..., False, False, False])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-ddc0a094-953d-4ee7-8741-9cf06bbd60b7' class='xr-section-summary-in' type='checkbox'  ><label for='section-ddc0a094-953d-4ee7-8741-9cf06bbd60b7' class='xr-section-summary' >Indexes: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-b9e29374-d3f0-4620-b9b0-32d2fbcc4c39' class='xr-index-data-in' type='checkbox'/><label for='index-b9e29374-d3f0-4620-b9b0-32d2fbcc4c39' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2015-05-06 00:00:00&#x27;, &#x27;2015-05-06 00:01:00&#x27;,
               &#x27;2015-05-06 00:02:00&#x27;, &#x27;2015-05-06 00:03:00&#x27;,
               &#x27;2015-05-06 00:04:00&#x27;, &#x27;2015-05-06 00:05:00&#x27;,
               &#x27;2015-05-06 00:06:00&#x27;, &#x27;2015-05-06 00:07:00&#x27;,
               &#x27;2015-05-06 00:08:00&#x27;, &#x27;2015-05-06 00:09:00&#x27;,
               ...
               &#x27;2015-05-11 23:50:00&#x27;, &#x27;2015-05-11 23:51:00&#x27;,
               &#x27;2015-05-11 23:52:00&#x27;, &#x27;2015-05-11 23:53:00&#x27;,
               &#x27;2015-05-11 23:54:00&#x27;, &#x27;2015-05-11 23:55:00&#x27;,
               &#x27;2015-05-11 23:56:00&#x27;, &#x27;2015-05-11 23:57:00&#x27;,
               &#x27;2015-05-11 23:58:00&#x27;, &#x27;2015-05-11 23:59:00&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=8640, freq=&#x27;T&#x27;))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-cc0855ba-f911-4a17-8b3e-082982ad4e8c' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-cc0855ba-f911-4a17-8b3e-082982ad4e8c' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>




```python
plt.plot(mordor.time,mordor.TP2_Wm2)
plt.plot(mordor.time[mordor.cs_mask],mordor.TP2_Wm2[mordor.cs_mask],ls='',marker='.')
```




    [<matplotlib.lines.Line2D at 0x7ff685cadc90>]




    
![png](calibration_melcol_cc_output_1.png)
    


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
Dimensions:          (station: 100, maintenancetime: 50, time: 518400)
Coordinates:
  * station          (station) int64 1 2 3 4 5 6 7 8 ... 94 95 96 97 98 99 100
  * maintenancetime  (maintenancetime) datetime64[ns] 2015-05-12T07:55:50 ......
  * time             (time) datetime64[ns] 2015-05-06 ... 2015-05-11T23:59:59
Data variables:
    ghi              (time, station) float64 0.0 nan nan 0.0 ... nan 0.0 0.0 nan
    szen             (time, station) float64 111.0 nan nan ... 109.4 109.4 nan
Attributes: (12/31)
    Conventions:               CF-1.10, ACDD-1.3
    title:                     TROPOS pyranometer network (PyrNet) observatio...
    history:                   2024-11-04T23:59:36: Merged level l1b by pyrne...
    institution:               Leibniz Institute for Tropospheric Research (T...
    source:                    TROPOS pyranometer network (PyrNet)
    references:                https://doi.org/10.5194/amt-9-1153-2016
    ...                        ...
    geospatial_lon_units:      degE
    time_coverage_start:       2015-05-06T00:00:00
    time_coverage_end:         2015-05-06T23:59:59
    time_coverage_duration:    P0DT23H59M59S
    time_coverage_resolution:  P0DT0H0M1S
    site:                      [&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;...</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-b849c982-d2b5-409b-8690-3beb7d39a334' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-b849c982-d2b5-409b-8690-3beb7d39a334' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>station</span>: 100</li><li><span class='xr-has-index'>maintenancetime</span>: 50</li><li><span class='xr-has-index'>time</span>: 518400</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-54abbe3f-c3fd-445b-b56f-92c376b7259e' class='xr-section-summary-in' type='checkbox'  checked><label for='section-54abbe3f-c3fd-445b-b56f-92c376b7259e' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>station</span></div><div class='xr-var-dims'>(station)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 2 3 4 5 6 ... 95 96 97 98 99 100</div><input id='attrs-697fc3f1-86e7-481f-b24e-30614532d5b1' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-697fc3f1-86e7-481f-b24e-30614532d5b1' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-e23ef5e7-291f-4fc6-991a-6d3906e9da86' class='xr-var-data-in' type='checkbox'><label for='data-e23ef5e7-291f-4fc6-991a-6d3906e9da86' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
        15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
        29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
        43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
        57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
        71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
        85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
        99, 100])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>maintenancetime</span></div><div class='xr-var-dims'>(maintenancetime)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-12T07:55:50 ... 2015-05-...</div><input id='attrs-9d5873ad-2633-40e5-968f-2ef98fae9362' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-9d5873ad-2633-40e5-968f-2ef98fae9362' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-479ae9e5-5754-4a59-bd4f-4b9d2b7f9f14' class='xr-var-data-in' type='checkbox'><label for='data-479ae9e5-5754-4a59-bd4f-4b9d2b7f9f14' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-12T07:55:50.000000000&#x27;, &#x27;2015-05-12T07:56:52.000000000&#x27;,
       &#x27;2015-05-12T07:57:20.000000000&#x27;, &#x27;2015-05-12T07:57:49.000000000&#x27;,
       &#x27;2015-05-12T07:58:21.000000000&#x27;, &#x27;2015-05-12T07:59:12.000000000&#x27;,
       &#x27;2015-05-12T07:59:46.000000000&#x27;, &#x27;2015-05-12T08:00:18.000000000&#x27;,
       &#x27;2015-05-12T08:00:55.000000000&#x27;, &#x27;2015-05-12T08:01:27.000000000&#x27;,
       &#x27;2015-05-12T08:02:19.000000000&#x27;, &#x27;2015-05-12T08:02:24.000000000&#x27;,
       &#x27;2015-05-12T08:02:59.000000000&#x27;, &#x27;2015-05-12T08:03:40.000000000&#x27;,
       &#x27;2015-05-12T08:04:14.000000000&#x27;, &#x27;2015-05-12T08:04:40.000000000&#x27;,
       &#x27;2015-05-12T08:05:06.000000000&#x27;, &#x27;2015-05-12T08:05:39.000000000&#x27;,
       &#x27;2015-05-12T08:06:44.000000000&#x27;, &#x27;2015-05-12T08:07:41.000000000&#x27;,
       &#x27;2015-05-12T08:08:21.000000000&#x27;, &#x27;2015-05-12T08:08:25.000000000&#x27;,
       &#x27;2015-05-12T08:08:44.000000000&#x27;, &#x27;2015-05-12T08:09:06.000000000&#x27;,
       &#x27;2015-05-12T08:09:08.000000000&#x27;, &#x27;2015-05-12T08:09:26.000000000&#x27;,
       &#x27;2015-05-12T08:09:35.000000000&#x27;, &#x27;2015-05-12T08:09:45.000000000&#x27;,
       &#x27;2015-05-12T08:10:08.000000000&#x27;, &#x27;2015-05-12T08:10:29.000000000&#x27;,
       &#x27;2015-05-12T08:10:52.000000000&#x27;, &#x27;2015-05-12T08:11:06.000000000&#x27;,
       &#x27;2015-05-12T08:11:17.000000000&#x27;, &#x27;2015-05-12T08:11:40.000000000&#x27;,
       &#x27;2015-05-12T08:11:42.000000000&#x27;, &#x27;2015-05-12T08:12:07.000000000&#x27;,
       &#x27;2015-05-12T08:12:14.000000000&#x27;, &#x27;2015-05-12T08:12:16.000000000&#x27;,
       &#x27;2015-05-12T08:12:39.000000000&#x27;, &#x27;2015-05-12T08:12:46.000000000&#x27;,
       &#x27;2015-05-12T08:12:59.000000000&#x27;, &#x27;2015-05-12T08:13:13.000000000&#x27;,
       &#x27;2015-05-12T08:13:17.000000000&#x27;, &#x27;2015-05-12T08:13:34.000000000&#x27;,
       &#x27;2015-05-12T08:13:55.000000000&#x27;, &#x27;2015-05-12T08:14:13.000000000&#x27;,
       &#x27;2015-05-12T08:14:34.000000000&#x27;, &#x27;2015-05-12T09:18:31.000000000&#x27;,
       &#x27;2015-05-12T09:19:12.000000000&#x27;, &#x27;2015-05-12T09:19:16.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-06 ... 2015-05-11T23:59:59</div><input id='attrs-4162c6cf-4a73-4c24-84d5-35e0af20b776' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-4162c6cf-4a73-4c24-84d5-35e0af20b776' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-72d98608-b6af-482c-9ee3-14c210f87b90' class='xr-var-data-in' type='checkbox'><label for='data-72d98608-b6af-482c-9ee3-14c210f87b90' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-06T00:00:00.000000000&#x27;, &#x27;2015-05-06T00:00:01.000000000&#x27;,
       &#x27;2015-05-06T00:00:02.000000000&#x27;, ..., &#x27;2015-05-11T23:59:57.000000000&#x27;,
       &#x27;2015-05-11T23:59:58.000000000&#x27;, &#x27;2015-05-11T23:59:59.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-325dd668-3259-4606-843e-6d5dc262588c' class='xr-section-summary-in' type='checkbox'  checked><label for='section-325dd668-3259-4606-843e-6d5dc262588c' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>ghi</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0 nan nan 0.0 ... nan 0.0 0.0 nan</div><input id='attrs-4dfa7f54-21fd-41dd-ac2b-53fb49b0d85a' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-4dfa7f54-21fd-41dd-ac2b-53fb49b0d85a' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-895041f8-45cd-4d61-b3e0-af2278b830db' class='xr-var-data-in' type='checkbox'><label for='data-895041f8-45cd-4d61-b3e0-af2278b830db' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>W m-2</dd><dt><span>long_name :</span></dt><dd>downwelling shortwave flux</dd><dt><span>standard_name :</span></dt><dd>downwelling_shortwave_flux_in_air</dd><dt><span>valid_range :</span></dt><dd>[    0 60000]</dd><dt><span>ancillary_variables :</span></dt><dd>maintenance_flag_ghi qc_flag_ghi</dd><dt><span>serial :</span></dt><dd>[&#x27;S12128.001&#x27;, &#x27;S12128.004&#x27;, &#x27;S12128.005&#x27;, &#x27;S12078.061&#x27;, &#x27;S12128.015&#x27;, &#x27;S12128.027&#x27;, &#x27;S12128.019&#x27;, &#x27;S12128.021&#x27;, &#x27;S12128.022&#x27;, &#x27;S12128.026&#x27;, &#x27;S12128.029&#x27;, &#x27;S12078.060&#x27;, &#x27;S12128.030&#x27;, &#x27;S12128.034&#x27;, &#x27;S12128.035&#x27;, &#x27;S12128.037&#x27;, &#x27;S12128.040&#x27;, &#x27;S12128.042&#x27;, &#x27;S12128.043&#x27;, &#x27;S12128.046&#x27;, &#x27;S12128.049&#x27;, &#x27;S12128.050&#x27;, &#x27;S12137.001&#x27;, &#x27;S12137.003&#x27;, &#x27;S12137.004&#x27;, &#x27;S12137.005&#x27;, &#x27;S12137.007&#x27;, &#x27;S12137.012&#x27;, &#x27;S12137.013&#x27;, &#x27;S12137.014&#x27;, &#x27;S12137.018&#x27;, &#x27;S12137.021&#x27;, &#x27;S12137.022&#x27;, &#x27;S12137.024&#x27;, &#x27;S12137.025&#x27;, &#x27;S12137.027&#x27;, &#x27;S12137.028&#x27;, &#x27;S12137.030&#x27;, &#x27;S12137.031&#x27;, &#x27;S12137.034&#x27;, &#x27;S12137.037&#x27;, &#x27;S12137.038&#x27;, &#x27;S12137.039&#x27;, &#x27;S12137.040&#x27;, &#x27;S12137.041&#x27;, &#x27;S12137.042&#x27;, &#x27;S12137.044&#x27;, &#x27;S12137.045&#x27;, &#x27;S12137.046&#x27;, &#x27;S12137.048&#x27;, &#x27;S12137.049&#x27;]</dd><dt><span>calibration_Cabsolute :</span></dt><dd>[142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
 142857.14285714 142857.14285714 142857.14285714 142857.14285714
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
 0 0 0 0 0 0 0 0 0 0 0 0 0 0]</dd><dt><span>hangle :</span></dt><dd>[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 0 0]</dd></dl></div><div class='xr-var-data'><pre>array([[ 0., nan, nan, ..., nan,  0., nan],
       [ 0., nan, nan, ..., nan,  0., nan],
       [ 0., nan, nan, ..., nan,  0., nan],
       ...,
       [ 0., nan, nan, ...,  0.,  0., nan],
       [ 0., nan, nan, ...,  0.,  0., nan],
       [ 0., nan, nan, ...,  0.,  0., nan]])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>szen</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>111.0 nan nan ... 109.4 109.4 nan</div><input id='attrs-7d6bfef3-93f9-4e97-bd47-ed4876ed38e0' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-7d6bfef3-93f9-4e97-bd47-ed4876ed38e0' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b9c6d05f-2b8e-4a59-a2c6-bdfef333cbf2' class='xr-var-data-in' type='checkbox'><label for='data-b9c6d05f-2b8e-4a59-a2c6-bdfef333cbf2' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>standard_name :</span></dt><dd>solar_zenith_angle</dd><dt><span>units :</span></dt><dd>degree</dd><dt><span>valid_range :</span></dt><dd>[    0 36000]</dd></dl></div><div class='xr-var-data'><pre>array([[111.03499603,          nan,          nan, ...,          nan,
        111.03499603,          nan],
       [111.03499603,          nan,          nan, ...,          nan,
        111.03499603,          nan],
       [111.03499603,          nan,          nan, ...,          nan,
        111.03499603,          nan],
       ...,
       [109.43000031,          nan,          nan, ..., 109.43000031,
        109.43000031,          nan],
       [109.43000031,          nan,          nan, ..., 109.43000031,
        109.43000031,          nan],
       [109.43000031,          nan,          nan, ..., 109.43000031,
        109.43000031,          nan]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-77c86e3e-13eb-4398-8d19-a597887bd939' class='xr-section-summary-in' type='checkbox'  ><label for='section-77c86e3e-13eb-4398-8d19-a597887bd939' class='xr-section-summary' >Indexes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>station</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-a5959465-fc63-488c-8fe5-a010c24f1fac' class='xr-index-data-in' type='checkbox'/><label for='index-a5959465-fc63-488c-8fe5-a010c24f1fac' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
        15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
        29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
        43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
        57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
        71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
        85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
        99, 100],
      dtype=&#x27;int64&#x27;, name=&#x27;station&#x27;))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>maintenancetime</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-16c5f363-e0f8-4887-b2ea-2aad4b89d1fb' class='xr-index-data-in' type='checkbox'/><label for='index-16c5f363-e0f8-4887-b2ea-2aad4b89d1fb' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2015-05-12 07:55:50&#x27;, &#x27;2015-05-12 07:56:52&#x27;,
               &#x27;2015-05-12 07:57:20&#x27;, &#x27;2015-05-12 07:57:49&#x27;,
               &#x27;2015-05-12 07:58:21&#x27;, &#x27;2015-05-12 07:59:12&#x27;,
               &#x27;2015-05-12 07:59:46&#x27;, &#x27;2015-05-12 08:00:18&#x27;,
               &#x27;2015-05-12 08:00:55&#x27;, &#x27;2015-05-12 08:01:27&#x27;,
               &#x27;2015-05-12 08:02:19&#x27;, &#x27;2015-05-12 08:02:24&#x27;,
               &#x27;2015-05-12 08:02:59&#x27;, &#x27;2015-05-12 08:03:40&#x27;,
               &#x27;2015-05-12 08:04:14&#x27;, &#x27;2015-05-12 08:04:40&#x27;,
               &#x27;2015-05-12 08:05:06&#x27;, &#x27;2015-05-12 08:05:39&#x27;,
               &#x27;2015-05-12 08:06:44&#x27;, &#x27;2015-05-12 08:07:41&#x27;,
               &#x27;2015-05-12 08:08:21&#x27;, &#x27;2015-05-12 08:08:25&#x27;,
               &#x27;2015-05-12 08:08:44&#x27;, &#x27;2015-05-12 08:09:06&#x27;,
               &#x27;2015-05-12 08:09:08&#x27;, &#x27;2015-05-12 08:09:26&#x27;,
               &#x27;2015-05-12 08:09:35&#x27;, &#x27;2015-05-12 08:09:45&#x27;,
               &#x27;2015-05-12 08:10:08&#x27;, &#x27;2015-05-12 08:10:29&#x27;,
               &#x27;2015-05-12 08:10:52&#x27;, &#x27;2015-05-12 08:11:06&#x27;,
               &#x27;2015-05-12 08:11:17&#x27;, &#x27;2015-05-12 08:11:40&#x27;,
               &#x27;2015-05-12 08:11:42&#x27;, &#x27;2015-05-12 08:12:07&#x27;,
               &#x27;2015-05-12 08:12:14&#x27;, &#x27;2015-05-12 08:12:16&#x27;,
               &#x27;2015-05-12 08:12:39&#x27;, &#x27;2015-05-12 08:12:46&#x27;,
               &#x27;2015-05-12 08:12:59&#x27;, &#x27;2015-05-12 08:13:13&#x27;,
               &#x27;2015-05-12 08:13:17&#x27;, &#x27;2015-05-12 08:13:34&#x27;,
               &#x27;2015-05-12 08:13:55&#x27;, &#x27;2015-05-12 08:14:13&#x27;,
               &#x27;2015-05-12 08:14:34&#x27;, &#x27;2015-05-12 09:18:31&#x27;,
               &#x27;2015-05-12 09:19:12&#x27;, &#x27;2015-05-12 09:19:16&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;maintenancetime&#x27;, freq=None))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-45c0c6b3-6fb0-40c9-8b4c-a8220fe2a167' class='xr-index-data-in' type='checkbox'/><label for='index-45c0c6b3-6fb0-40c9-8b4c-a8220fe2a167' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2015-05-06 00:00:00&#x27;, &#x27;2015-05-06 00:00:01&#x27;,
               &#x27;2015-05-06 00:00:02&#x27;, &#x27;2015-05-06 00:00:03&#x27;,
               &#x27;2015-05-06 00:00:04&#x27;, &#x27;2015-05-06 00:00:05&#x27;,
               &#x27;2015-05-06 00:00:06&#x27;, &#x27;2015-05-06 00:00:07&#x27;,
               &#x27;2015-05-06 00:00:08&#x27;, &#x27;2015-05-06 00:00:09&#x27;,
               ...
               &#x27;2015-05-11 23:59:50&#x27;, &#x27;2015-05-11 23:59:51&#x27;,
               &#x27;2015-05-11 23:59:52&#x27;, &#x27;2015-05-11 23:59:53&#x27;,
               &#x27;2015-05-11 23:59:54&#x27;, &#x27;2015-05-11 23:59:55&#x27;,
               &#x27;2015-05-11 23:59:56&#x27;, &#x27;2015-05-11 23:59:57&#x27;,
               &#x27;2015-05-11 23:59:58&#x27;, &#x27;2015-05-11 23:59:59&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=518400, freq=None))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-0087d5a5-24f0-4976-8e29-94c49788815c' class='xr-section-summary-in' type='checkbox'  ><label for='section-0087d5a5-24f0-4976-8e29-94c49788815c' class='xr-section-summary' >Attributes: <span>(31)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>Conventions :</span></dt><dd>CF-1.10, ACDD-1.3</dd><dt><span>title :</span></dt><dd>TROPOS pyranometer network (PyrNet) observational data set</dd><dt><span>history :</span></dt><dd>2024-11-04T23:59:36: Merged level l1b by pyrnet version 0.2.16; </dd><dt><span>institution :</span></dt><dd>Leibniz Institute for Tropospheric Research (TROPOS)</dd><dt><span>source :</span></dt><dd>TROPOS pyranometer network (PyrNet)</dd><dt><span>references :</span></dt><dd>https://doi.org/10.5194/amt-9-1153-2016</dd><dt><span>Department :</span></dt><dd>Remote Sensing of Atmospheric Processes</dd><dt><span>Department_team :</span></dt><dd>Clouds, Aerosol and Radiation</dd><dt><span>Address :</span></dt><dd>Permoser Str. 15, 04318 Leipzig, Germany</dd><dt><span>Contact_person :</span></dt><dd>Andreas Macke and the clouds, aerosol and radiation team of the remote sensing department, mailto:andreas.macke@tropos.de</dd><dt><span>Contributor_name :</span></dt><dd>[&#x27;Leibniz-Institute for Tropospheric Research (TROPOS)&#x27;, &#x27;Andreas Macke (TROPOS)&#x27;, &#x27;Hartwig Deneke (TROPOS)&#x27;, &#x27;Bomidi Lakshmi Madhavan (TROPOS)&#x27;, &#x27;Jonas Witthuhn (TROPOS)&#x27;, &#x27;Ronny Badeke (TROPOS)&#x27;]</dd><dt><span>Contributor_role :</span></dt><dd>[&#x27;supervision, setup, maintenance, tear-down, data-processing, site-planning&#x27;, &#x27;project-coordinaton, supervision,&#x27;, &#x27;supervision&#x27;, &#x27;supervision, data-processing, quality-control, site-planning&#x27;, &#x27;supervision, data-processing, quality-control, setup, tear-down, maintenance&#x27;, &#x27;setup, tear-down, maintenance&#x27;]</dd><dt><span>Authors_software :</span></dt><dd>Hartwig Deneke, Jonas Witthuhn, mailto:deneke@tropos.de</dd><dt><span>Creator_name :</span></dt><dd>Jonas Witthuhn</dd><dt><span>Project :</span></dt><dd>Intensive Measurement Campaign Melpitz Column (MelCol)</dd><dt><span>Standard_name_vocabulary :</span></dt><dd>CF Standard Name Table v81</dd><dt><span>License :</span></dt><dd>CC-BY-SA 3.0</dd><dt><span>processing_level :</span></dt><dd>l1b</dd><dt><span>product_version :</span></dt><dd>0.2.16</dd><dt><span>date_created :</span></dt><dd>2024-11-05T00:01:39</dd><dt><span>geospatial_lat_min :</span></dt><dd>51.5250188772739</dd><dt><span>geospatial_lat_max :</span></dt><dd>51.5253452457728</dd><dt><span>geospatial_lat_units :</span></dt><dd>degN</dd><dt><span>geospatial_lon_min :</span></dt><dd>12.926223022881615</dd><dt><span>geospatial_lon_max :</span></dt><dd>12.926649323253269</dd><dt><span>geospatial_lon_units :</span></dt><dd>degE</dd><dt><span>time_coverage_start :</span></dt><dd>2015-05-06T00:00:00</dd><dt><span>time_coverage_end :</span></dt><dd>2015-05-06T23:59:59</dd><dt><span>time_coverage_duration :</span></dt><dd>P0DT23H59M59S</dd><dt><span>time_coverage_resolution :</span></dt><dd>P0DT0H0M1S</dd><dt><span>site :</span></dt><dd>[&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;]</dd></dl></div></li></ul></div></div>



### Calibration
The calibration follows the [ISO 9847:1992 - Solar energy — Calibration of field pyranometers by comparison to a reference pyranometer](https://archive.org/details/gov.in.is.iso.9847.1992). For selecting and masking the data.
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
# Set flux values to nan if no pyranometer is installed.
pyr.ghi.values = pyr.ghi.where(mainmask).values

# convert to measured voltage
pyr.ghi.values = pyr.ghi.values * 7 * 1e-6

# Step 1, select data
pyr = pyr.where(pyr.szen<80, drop=True)
pyr.ghi.values = pyr.ghi.where(pyr.ghi>0.033/300.).values

```

#### Step 2
Apply new absolute calibration from `calibration_melcol.ipynb`


```python
calib_new = pyrnet.utils.read_json("pyrnet_calib_new.json")
calib_new = calib_new[list(calib_new.keys())[0]]
for station in pyr.station:
    pyr.ghi.sel(station=station).values /= calib_new[f"{station:03d}"][0] * 1e-6 


```

#### Step 3
Interpolate reference to PyrNet samples and combine to a single Dataset


```python
# interpolate reference to PyrNet
mordor = mordor.interp(time=pyr.time).interpolate_na()


# Calibration datasets for main and extra pyranometer
Cds_main = xr.Dataset(
    data_vars={
        'reference_Wm2': ('time', mordor['TP2_Wm2'].data[mordor.cs_mask]),
        'pyrnet_Wm2': (('time','station'), pyr['ghi'].data[mordor.cs_mask,:]),
        'szen': ('time',pyr.szen.mean("station",skipna=True).data[mordor.cs_mask])
    },
    coords= {
        "time": pyr.time[mordor.cs_mask],
        "station": pyr.station
    }
)

```

#### Step 4
Remove outliers from series using xarray grouping and apply function.


```python
def remove_outliers(x):
    """
    x is an xarray dataset containing these variables:
    coords: 'time' - datetime64
    'pyrnet_V' - array - voltage measures of pyranometer
    'reference_Wm2' - array - measured irradiance of reference
    """

    # calculate calibration series for single samples
    C = x['pyrnet_Wm2'] / x['reference_Wm2']
    # integrated series 
    ix = x.integrate('time')
    M = ix['pyrnet_Wm2'] / ix['reference_Wm2']
    
    while np.any(np.abs(C-M) > 0.01*M):
        #calculate as long there are outliers deviating more than 2 percent
        x = x.where(np.abs(C-M) < 0.01*M)
        C = x['pyrnet_Wm2'] / x['reference_Wm2']
        #integrated series 
        ix = x.integrate('time')
        M = ix['pyrnet_Wm2'] / ix['reference_Wm2']
        
    #return the reduced dataset x
    return x

# remove outliers
Cds_main = Cds_main.groupby('time.hour').apply(remove_outliers)

# hourly mean
Cds_main = Cds_main.coarsen(time=60*60,boundary='trim').mean(skipna=True)

```

#### Step 5
The series of measured reference and pyrnet irradiance is now without outliers. Now we can fit a cubic function - (reference/pyrnet) vs. szen - to determine the cosine correction function for pyrnet.


```python
fig,ax = plt.subplots(1,1)
p = ax.plot(Cds_main.szen,Cds_main.reference_Wm2/Cds_main.pyrnet_Wm2,ls='',marker='.')
```


    
![png](calibration_melcol_cc_output_0.png)
    



```python
def cubic_function(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d

def objective_function(coefficients, x, y):
    a, b, c, d = coefficients
    y_pred = cubic_function(x, a, b, c, d)
    error = np.sum((y - y_pred) ** 2)
    return error
bounds = [(-5, 5), (-5, 5), (-5, 5), (-5, 5)]


ratio = Cds_main.reference_Wm2/Cds_main.pyrnet_Wm2

result = differential_evolution(
    objective_function,
    args=(np.cos(np.deg2rad(Cds_main.szen)), ratio),
    bounds=bounds,
    seed=1
)
result
```




     message: Optimization terminated successfully.
     success: True
         fun: 0.39807904993609633
           x: [-2.335e+00  4.402e+00 -2.437e+00  1.379e+00]
         nit: 56
        nfev: 3545
         jac: [ 2.380e-04  4.970e-04  7.820e-04  1.921e-04]



### Results


```python
print("Best cubic fit:")
a3, a2, a1, a0 = result.x
display(Latex(
    rf"""
    {a3:+.3f}x^3{a2:+.3f}x^2{a1:+.3f}x{a0:+.3f}
    """
))
```

    Best cubic fit:




    -2.335x^3+4.402x^2-2.437x+1.379
    



```python
# Coefficients from other calibrations:

display(Latex(
    rf"""
    MelCol: {a3:+.3f}x^3{a2:+.3f}x^2{a1:+.3f}x{a0:+.3f}
    """
))

b3, b2, b1, b0 = [ -3.01, 5.59, -3.04, 1.45 ]
display(Latex(
    rf"""
    MetPVNet: {b3:+.3f}x^3{b2:+.3f}x^2{b1:+.3f}x{b0:+.3f}
    """
))
c3, c2, c1, c0 = [ -2.227, 4.366, -2.524, 1.385 ]
display(Latex(
    rf"""
    S2VSR: {c3:+.3f}x^3{c2:+.3f}x^2{c1:+.3f}x{c0:+.3f}
    """
))
```



    MelCol: -2.335x^3+4.402x^2-2.437x+1.379
    




    MetPVNet: -3.010x^3+5.590x^2-3.040x+1.450
    




    S2VSR: -2.227x^3+4.366x^2-2.524x+1.385
    



```python
szen = np.arange(1,80)
mu0 = np.cos(np.deg2rad(szen))

fig,ax = plt.subplots(1,1)
p = ax.plot(Cds_main.szen,Cds_main.reference_Wm2/Cds_main.pyrnet_Wm2,ls='',marker='.')
ax.plot(szen, a3*mu0**3 + a2*mu0**2 + a1*mu0 + a0,
         color='C1', label=f'best-cubic-fit: {a3:+.3f}x^3{a2:+.3f}x^2{a1:+.3f}x{a0:+.3f}')
ax.plot(szen, b3*mu0**3 + b2*mu0**2 + b1*mu0 + b0,
         color='C2', label=f'MetPvNet: {b3:+.3f}x^3{b2:+.3f}x^2{b1:+.3f}x{b0:+.3f}')
ax.plot(szen, c3*mu0**3 + c2*mu0**2 + c1*mu0 + c0,
         color='C3', label=f'S2VSR: {c3:+.3f}x^3{c2:+.3f}x^2{c1:+.3f}x{c0:+.3f}')

ax.legend()
ax.set_xlabel('Zenith angle [deg] ')
ax.set_ylabel('reference_ghi/pyr_ghi // correction_factor [-]')
ax.grid(True)


```


    
![png](calibration_melcol_cc_output_26_0.png)
    



```python
# dump to json
calibjson = {"2015-05-06": {"CC":list(result.x[::-1])}}
with open("pyrnet_calib_cc_new.json","w") as txt:
    json.dump(calibjson, txt)
```
