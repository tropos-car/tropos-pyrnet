## 2015 MelCol
```
processed with pyrnet-0.2.16
```

The PyrNet was setup for calibration in a dense array on the Melpitz measurement field from 2015-05-06 to 2015-05-11. Cross-calibration is done versus reference observations from the TROPOS MObile RaDiation ObseRvatory (MORDOR) station.

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
pf_mordor = "mordor/{date:%Y-%m-%d}_Radiation.dat"
pf_pyrnet = "l1b_network/{date:%Y-%m-%d}_P1D_pyrnet_melcol_n000l1bf1s.c01.nc"
dates = pd.date_range("2015-05-06","2015-05-11")
stations = np.arange(1,101)

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
    keep_vars = ['TP2_Wm2'] # global shortwave irradiance
    drop_vars = [v for v in dst if v not in keep_vars]
    dst = dst.drop_vars(drop_vars)
    
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
Dimensions:  (time: 5183535)
Coordinates:
  * time     (time) datetime64[ns] 2015-05-06 ... 2015-05-11T23:59:59.900000
Data variables:
    TP2_Wm2  (time) float64 0.0 0.0 0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0 0.0 0.0</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-63e440f1-6a9c-4372-b2f4-c91d7060fc11' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-63e440f1-6a9c-4372-b2f4-c91d7060fc11' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>time</span>: 5183535</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-bb59f3df-1f8d-4753-84ee-89d59723313b' class='xr-section-summary-in' type='checkbox'  checked><label for='section-bb59f3df-1f8d-4753-84ee-89d59723313b' class='xr-section-summary' >Coordinates: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-06 ... 2015-05-11T23:59:...</div><input id='attrs-b6416221-9124-457f-8f45-559e8b9c0851' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-b6416221-9124-457f-8f45-559e8b9c0851' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-dfacedaf-d55f-41d1-bf13-880c3b9c20c6' class='xr-var-data-in' type='checkbox'><label for='data-dfacedaf-d55f-41d1-bf13-880c3b9c20c6' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-06T00:00:00.000000000&#x27;, &#x27;2015-05-06T00:00:00.300000000&#x27;,
       &#x27;2015-05-06T00:00:00.400000000&#x27;, ..., &#x27;2015-05-11T23:59:59.700000000&#x27;,
       &#x27;2015-05-11T23:59:59.800000000&#x27;, &#x27;2015-05-11T23:59:59.900000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-025117ff-72eb-4a2e-ac2f-dcd58c2139df' class='xr-section-summary-in' type='checkbox'  checked><label for='section-025117ff-72eb-4a2e-ac2f-dcd58c2139df' class='xr-section-summary' >Data variables: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>TP2_Wm2</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0</div><input id='attrs-5cb69921-751c-464e-a8d4-8a9ebf8cc4b7' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-5cb69921-751c-464e-a8d4-8a9ebf8cc4b7' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-7179ab6d-0d8d-4fca-906c-c4039b723988' class='xr-var-data-in' type='checkbox'><label for='data-7179ab6d-0d8d-4fca-906c-c4039b723988' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([0., 0., 0., ..., 0., 0., 0.])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-04a92278-9972-4ad3-82d2-47f0582dd4d7' class='xr-section-summary-in' type='checkbox'  ><label for='section-04a92278-9972-4ad3-82d2-47f0582dd4d7' class='xr-section-summary' >Indexes: <span>(1)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-71c5678b-1b6d-4242-ae35-ded6efae872d' class='xr-index-data-in' type='checkbox'/><label for='index-71c5678b-1b6d-4242-ae35-ded6efae872d' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([       &#x27;2015-05-06 00:00:00&#x27;, &#x27;2015-05-06 00:00:00.300000&#x27;,
               &#x27;2015-05-06 00:00:00.400000&#x27;, &#x27;2015-05-06 00:00:00.500000&#x27;,
               &#x27;2015-05-06 00:00:00.600000&#x27;, &#x27;2015-05-06 00:00:00.700000&#x27;,
               &#x27;2015-05-06 00:00:00.800000&#x27;, &#x27;2015-05-06 00:00:00.900000&#x27;,
                      &#x27;2015-05-06 00:00:01&#x27;, &#x27;2015-05-06 00:00:01.100000&#x27;,
               ...
                      &#x27;2015-05-11 23:59:59&#x27;, &#x27;2015-05-11 23:59:59.100000&#x27;,
               &#x27;2015-05-11 23:59:59.200000&#x27;, &#x27;2015-05-11 23:59:59.300000&#x27;,
               &#x27;2015-05-11 23:59:59.400000&#x27;, &#x27;2015-05-11 23:59:59.500000&#x27;,
               &#x27;2015-05-11 23:59:59.600000&#x27;, &#x27;2015-05-11 23:59:59.700000&#x27;,
               &#x27;2015-05-11 23:59:59.800000&#x27;, &#x27;2015-05-11 23:59:59.900000&#x27;],
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=5183535, freq=None))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-1610aa68-82b8-46d4-a459-b4e0f4b4ca2b' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-1610aa68-82b8-46d4-a459-b4e0f4b4ca2b' class='xr-section-summary'  title='Expand/collapse section'>Attributes: <span>(0)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'></dl></div></li></ul></div></div>



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
    site:                      [&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;...</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-40ad30f8-0003-4056-bacc-fa4fb8a262c9' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-40ad30f8-0003-4056-bacc-fa4fb8a262c9' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>station</span>: 100</li><li><span class='xr-has-index'>maintenancetime</span>: 50</li><li><span class='xr-has-index'>time</span>: 518400</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-83e450b1-ca15-40db-8208-bead2cce3ed6' class='xr-section-summary-in' type='checkbox'  checked><label for='section-83e450b1-ca15-40db-8208-bead2cce3ed6' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>station</span></div><div class='xr-var-dims'>(station)</div><div class='xr-var-dtype'>int64</div><div class='xr-var-preview xr-preview'>1 2 3 4 5 6 ... 95 96 97 98 99 100</div><input id='attrs-29edbdb1-96b7-4910-ae2d-4ce50c3df13e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-29edbdb1-96b7-4910-ae2d-4ce50c3df13e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-30b77282-47f2-4d9d-bdb0-9098f7332c06' class='xr-var-data-in' type='checkbox'><label for='data-30b77282-47f2-4d9d-bdb0-9098f7332c06' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
        15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
        29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
        43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
        57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
        71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
        85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
        99, 100])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>maintenancetime</span></div><div class='xr-var-dims'>(maintenancetime)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-12T07:55:50 ... 2015-05-...</div><input id='attrs-f25cd103-790e-43b6-8724-b92fe0de97e2' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-f25cd103-790e-43b6-8724-b92fe0de97e2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-301d8ff3-e4db-4f6b-afba-da36336dd011' class='xr-var-data-in' type='checkbox'><label for='data-301d8ff3-e4db-4f6b-afba-da36336dd011' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-12T07:55:50.000000000&#x27;, &#x27;2015-05-12T07:56:52.000000000&#x27;,
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
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>time</span></div><div class='xr-var-dims'>(time)</div><div class='xr-var-dtype'>datetime64[ns]</div><div class='xr-var-preview xr-preview'>2015-05-06 ... 2015-05-11T23:59:59</div><input id='attrs-855c635d-f281-4b10-9ea0-9f4a0477085e' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-855c635d-f281-4b10-9ea0-9f4a0477085e' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-31fb070c-d66d-49cc-8be2-b89813c62896' class='xr-var-data-in' type='checkbox'><label for='data-31fb070c-d66d-49cc-8be2-b89813c62896' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;2015-05-06T00:00:00.000000000&#x27;, &#x27;2015-05-06T00:00:01.000000000&#x27;,
       &#x27;2015-05-06T00:00:02.000000000&#x27;, ..., &#x27;2015-05-11T23:59:57.000000000&#x27;,
       &#x27;2015-05-11T23:59:58.000000000&#x27;, &#x27;2015-05-11T23:59:59.000000000&#x27;],
      dtype=&#x27;datetime64[ns]&#x27;)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-08589d3b-4bda-4585-9d4c-93bcf4ce86d8' class='xr-section-summary-in' type='checkbox'  checked><label for='section-08589d3b-4bda-4585-9d4c-93bcf4ce86d8' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>ghi</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>0.0 nan nan 0.0 ... nan 0.0 0.0 nan</div><input id='attrs-d353805d-328f-4468-b1e6-fa57a9e770ff' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-d353805d-328f-4468-b1e6-fa57a9e770ff' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-b890e4eb-852a-4529-ae71-7560657c21a0' class='xr-var-data-in' type='checkbox'><label for='data-b890e4eb-852a-4529-ae71-7560657c21a0' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>W m-2</dd><dt><span>long_name :</span></dt><dd>downwelling shortwave flux</dd><dt><span>standard_name :</span></dt><dd>downwelling_shortwave_flux_in_air</dd><dt><span>valid_range :</span></dt><dd>[    0 60000]</dd><dt><span>ancillary_variables :</span></dt><dd>maintenance_flag_ghi qc_flag_ghi</dd><dt><span>serial :</span></dt><dd>[&#x27;S12128.001&#x27;, &#x27;S12128.004&#x27;, &#x27;S12128.005&#x27;, &#x27;S12078.061&#x27;, &#x27;S12128.015&#x27;, &#x27;S12128.027&#x27;, &#x27;S12128.019&#x27;, &#x27;S12128.021&#x27;, &#x27;S12128.022&#x27;, &#x27;S12128.026&#x27;, &#x27;S12128.029&#x27;, &#x27;S12078.060&#x27;, &#x27;S12128.030&#x27;, &#x27;S12128.034&#x27;, &#x27;S12128.035&#x27;, &#x27;S12128.037&#x27;, &#x27;S12128.040&#x27;, &#x27;S12128.042&#x27;, &#x27;S12128.043&#x27;, &#x27;S12128.046&#x27;, &#x27;S12128.049&#x27;, &#x27;S12128.050&#x27;, &#x27;S12137.001&#x27;, &#x27;S12137.003&#x27;, &#x27;S12137.004&#x27;, &#x27;S12137.005&#x27;, &#x27;S12137.007&#x27;, &#x27;S12137.012&#x27;, &#x27;S12137.013&#x27;, &#x27;S12137.014&#x27;, &#x27;S12137.018&#x27;, &#x27;S12137.021&#x27;, &#x27;S12137.022&#x27;, &#x27;S12137.024&#x27;, &#x27;S12137.025&#x27;, &#x27;S12137.027&#x27;, &#x27;S12137.028&#x27;, &#x27;S12137.030&#x27;, &#x27;S12137.031&#x27;, &#x27;S12137.034&#x27;, &#x27;S12137.037&#x27;, &#x27;S12137.038&#x27;, &#x27;S12137.039&#x27;, &#x27;S12137.040&#x27;, &#x27;S12137.041&#x27;, &#x27;S12137.042&#x27;, &#x27;S12137.044&#x27;, &#x27;S12137.045&#x27;, &#x27;S12137.046&#x27;, &#x27;S12137.048&#x27;, &#x27;S12137.049&#x27;]</dd><dt><span>calibration_Cabsolute :</span></dt><dd>[142857.14285714 142857.14285714 142857.14285714 142857.14285714
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
       [ 0., nan, nan, ...,  0.,  0., nan]])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>szen</span></div><div class='xr-var-dims'>(time, station)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>111.0 nan nan ... 109.4 109.4 nan</div><input id='attrs-847f2fc6-0f31-493c-af6f-b6fe8cad020c' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-847f2fc6-0f31-493c-af6f-b6fe8cad020c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-044e830b-8015-44ce-8bdd-6706d589ba0e' class='xr-var-data-in' type='checkbox'><label for='data-044e830b-8015-44ce-8bdd-6706d589ba0e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>standard_name :</span></dt><dd>solar_zenith_angle</dd><dt><span>units :</span></dt><dd>degree</dd><dt><span>valid_range :</span></dt><dd>[    0 36000]</dd></dl></div><div class='xr-var-data'><pre>array([[111.03499603,          nan,          nan, ...,          nan,
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
        109.43000031,          nan]])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-fd7d9837-e49f-49bc-bdb5-700540c2479c' class='xr-section-summary-in' type='checkbox'  ><label for='section-fd7d9837-e49f-49bc-bdb5-700540c2479c' class='xr-section-summary' >Indexes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>station</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-3f189f90-c10a-41b6-8d8a-160b86c06603' class='xr-index-data-in' type='checkbox'/><label for='index-3f189f90-c10a-41b6-8d8a-160b86c06603' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
        15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
        29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
        43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
        57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
        71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,
        85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
        99, 100],
      dtype=&#x27;int64&#x27;, name=&#x27;station&#x27;))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>maintenancetime</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-735c25ae-a89d-41bf-bf9c-b6b4b83d51e8' class='xr-index-data-in' type='checkbox'/><label for='index-735c25ae-a89d-41bf-bf9c-b6b4b83d51e8' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2015-05-12 07:55:50&#x27;, &#x27;2015-05-12 07:56:52&#x27;,
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
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;maintenancetime&#x27;, freq=None))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>time</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-7e8225d3-c01e-43de-aa4f-543a57948408' class='xr-index-data-in' type='checkbox'/><label for='index-7e8225d3-c01e-43de-aa4f-543a57948408' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(DatetimeIndex([&#x27;2015-05-06 00:00:00&#x27;, &#x27;2015-05-06 00:00:01&#x27;,
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
              dtype=&#x27;datetime64[ns]&#x27;, name=&#x27;time&#x27;, length=518400, freq=None))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-b094a366-86f7-497c-85b7-b7d402d528b9' class='xr-section-summary-in' type='checkbox'  ><label for='section-b094a366-86f7-497c-85b7-b7d402d528b9' class='xr-section-summary' >Attributes: <span>(31)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>Conventions :</span></dt><dd>CF-1.10, ACDD-1.3</dd><dt><span>title :</span></dt><dd>TROPOS pyranometer network (PyrNet) observational data set</dd><dt><span>history :</span></dt><dd>2024-11-04T23:59:36: Merged level l1b by pyrnet version 0.2.16; </dd><dt><span>institution :</span></dt><dd>Leibniz Institute for Tropospheric Research (TROPOS)</dd><dt><span>source :</span></dt><dd>TROPOS pyranometer network (PyrNet)</dd><dt><span>references :</span></dt><dd>https://doi.org/10.5194/amt-9-1153-2016</dd><dt><span>Department :</span></dt><dd>Remote Sensing of Atmospheric Processes</dd><dt><span>Department_team :</span></dt><dd>Clouds, Aerosol and Radiation</dd><dt><span>Address :</span></dt><dd>Permoser Str. 15, 04318 Leipzig, Germany</dd><dt><span>Contact_person :</span></dt><dd>Andreas Macke and the clouds, aerosol and radiation team of the remote sensing department, mailto:andreas.macke@tropos.de</dd><dt><span>Contributor_name :</span></dt><dd>[&#x27;Leibniz-Institute for Tropospheric Research (TROPOS)&#x27;, &#x27;Andreas Macke (TROPOS)&#x27;, &#x27;Hartwig Deneke (TROPOS)&#x27;, &#x27;Bomidi Lakshmi Madhavan (TROPOS)&#x27;, &#x27;Jonas Witthuhn (TROPOS)&#x27;, &#x27;Ronny Badeke (TROPOS)&#x27;]</dd><dt><span>Contributor_role :</span></dt><dd>[&#x27;supervision, setup, maintenance, tear-down, data-processing, site-planning&#x27;, &#x27;project-coordinaton, supervision,&#x27;, &#x27;supervision&#x27;, &#x27;supervision, data-processing, quality-control, site-planning&#x27;, &#x27;supervision, data-processing, quality-control, setup, tear-down, maintenance&#x27;, &#x27;setup, tear-down, maintenance&#x27;]</dd><dt><span>Authors_software :</span></dt><dd>Hartwig Deneke, Jonas Witthuhn, mailto:deneke@tropos.de</dd><dt><span>Creator_name :</span></dt><dd>Jonas Witthuhn</dd><dt><span>Project :</span></dt><dd>Intensive Measurement Campaign Melpitz Column (MelCol)</dd><dt><span>Standard_name_vocabulary :</span></dt><dd>CF Standard Name Table v81</dd><dt><span>License :</span></dt><dd>CC-BY-SA 3.0</dd><dt><span>processing_level :</span></dt><dd>l1b</dd><dt><span>product_version :</span></dt><dd>0.2.16</dd><dt><span>date_created :</span></dt><dd>2024-11-05T00:01:39</dd><dt><span>geospatial_lat_min :</span></dt><dd>51.5250188772739</dd><dt><span>geospatial_lat_max :</span></dt><dd>51.5253452457728</dd><dt><span>geospatial_lat_units :</span></dt><dd>degN</dd><dt><span>geospatial_lon_min :</span></dt><dd>12.926223022881615</dd><dt><span>geospatial_lon_max :</span></dt><dd>12.926649323253269</dd><dt><span>geospatial_lon_units :</span></dt><dd>degE</dd><dt><span>time_coverage_start :</span></dt><dd>2015-05-06T00:00:00</dd><dt><span>time_coverage_end :</span></dt><dd>2015-05-06T23:59:59</dd><dt><span>time_coverage_duration :</span></dt><dd>P0DT23H59M59S</dd><dt><span>time_coverage_resolution :</span></dt><dd>P0DT0H0M1S</dd><dt><span>site :</span></dt><dd>[&#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;, &#x27;&#x27;]</dd></dl></div></li></ul></div></div>



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
# Set flux values to nan if no pyranometer is installed.
pyr.ghi.values = pyr.ghi.where(mainmask).values

# convert to measured voltage
pyr.ghi.values = pyr.ghi.values * 7 * 1e-6

# Step 1, select data
pyr = pyr.where(pyr.szen<80, drop=True)
pyr.ghi.values = pyr.ghi.where(pyr.ghi>0.033/300.).values

```

#### Step 2
Interpolate reference to PyrNet samples and combine to a single Dataset


```python
# interpolate reference to PyrNet
mordor = mordor.interp(time=pyr.time)

# Calibration datasets for main and extra pyranometer
Cds_main = xr.Dataset(
    data_vars={
        'reference_Wm2': ('time', mordor['TP2_Wm2'].data),
        'pyrnet_V': (('time','station'), pyr['ghi'].data)
    },
    coords= {
        "time": pyr.time,
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
Cds_main = Cds_main.coarsen(time=60*60,boundary='trim').mean(skipna=True)

```

#### Step 4
The series of measured voltage and irradiance is now without outliers. So we use equation 1 again to calculate from this reduced series the calibration factor for the instant samples.


```python
C_main = 1e6*Cds_main['pyrnet_V'] / Cds_main['reference_Wm2']
```

#### Step 5
We just found the Calibration factor to be the mean of the reduced calibration factor series and the uncertainty to be the standard deviation of this reduced series. Steo 3, 4 and 5 are done for every pyranometer seperate.


```python
C_main_mean = C_main.mean(dim='time',skipna=True)
C_main_std = C_main.std(dim='time',skipna=True)
```

### Results


```python
#|dropcode
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_title("Main Pyranometer")
ax.plot(C_main.time, C_main, ls ="", marker='.')
ax.set_xlabel("Date")
ax.set_ylabel("Calibration factor (uV / Wm-2)")
ax.grid(True)
fig.show()

plt.figure()
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_title("Main Pyranometer")
ax.plot(pyr.szen.interp_like(C_main), C_main, ls ="", marker='.')
ax.set_xlabel("solar zenith angle (deg)")
ax.set_ylabel("Calibration factor (uV / Wm-2)")
ax.grid(True)
fig.show()



```


    
![png](calibration_melcol_output_0.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](calibration_melcol_output_2.png)
    



```python
calibration_new = {}
print(f"Box:    Main       ,     Extra  ")
for box in C_main_mean.station:
    Cm = float(C_main_mean.sel(station=box).values)
    Um = float(C_main_std.sel(station=box).values)
    
    calibration_new.update({
        f"{box:03d}": [np.round(Cm,2), None]
    })
    print(f"{box:3d}: {Cm:.2f} +- {Um:.3f} , {None}")

calibjson = {"2015-05-06": calibration_new}
with open("pyrnet_calib_new.json","w") as txt:
    json.dump(calibjson, txt)
```

    Box:    Main       ,     Extra  
      1: 7.44 +- 0.204 , None
      4: 7.41 +- 0.210 , None
      5: 7.40 +- 0.182 , None
      6: 6.54 +- 0.088 , None
     15: 7.34 +- 0.181 , None
     16: 7.59 +- 0.190 , None
     19: 7.46 +- 0.262 , None
     21: 7.40 +- 0.200 , None
     22: 7.40 +- 0.205 , None
     26: 7.49 +- 0.150 , None
     28: 7.47 +- 0.244 , None
     29: 7.25 +- 0.165 , None
     30: 7.55 +- 0.176 , None
     34: 6.93 +- 0.162 , None
     35: 7.43 +- 0.182 , None
     37: 7.52 +- 0.148 , None
     40: 7.37 +- 0.199 , None
     42: 7.53 +- 0.213 , None
     43: 7.26 +- 0.170 , None
     46: 7.65 +- 0.179 , None
     49: 7.63 +- 0.241 , None
     50: 7.59 +- 0.176 , None
     51: 7.35 +- 0.087 , None
     53: 7.36 +- 0.256 , None
     54: 7.35 +- 0.197 , None
     55: 7.18 +- 0.153 , None
     57: 6.63 +- 0.170 , None
     62: 7.24 +- 0.167 , None
     63: 7.38 +- 0.118 , None
     64: 7.19 +- 0.111 , None
     68: 6.68 +- 0.082 , None
     71: 7.32 +- 0.186 , None
     72: 7.26 +- 0.184 , None
     74: 7.41 +- 0.096 , None
     75: 6.52 +- 0.194 , None
     77: 7.36 +- 0.089 , None
     78: 6.63 +- 0.154 , None
     80: 6.85 +- 0.102 , None
     81: 7.34 +- 0.088 , None
     84: 7.07 +- 0.198 , None
     87: 7.26 +- 0.136 , None
     88: 6.91 +- 0.102 , None
     89: 7.31 +- 0.132 , None
     90: 7.17 +- 0.152 , None
     91: 7.37 +- 0.182 , None
     92: 7.17 +- 0.169 , None
     94: 6.57 +- 0.136 , None
     95: 7.32 +- 0.203 , None
     96: 7.39 +- 0.149 , None
     98: 7.20 +- 0.210 , None
     99: 7.34 +- 0.126 , None

