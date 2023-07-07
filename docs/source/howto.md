---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


# How to ... 

## ... get data from tds.tropos.de?
Currently, the processed PyrNet data is hosted on the [TROPOS Thredds server](https://tds.tropos.de/thredds/catalog/scccJher/catalog.html).
Data can be obtained, e.g.,  via the OPeNDAP protocol, or using the [pyrnet python module](https://github.com/jonas-witthuhn/PyrNet). 

```{code-cell}
:tags: [hide-input]

import datetime as dt
from pyrnet import pyrnet

print(pyrnet.read_thredds.__doc__)
```

For example, read and merge two days from the S2VSR campaign:
```{code-cell}
:tags: [hide-input]

campaign='s2vsr'
dates = [dt.datetime(2023,6,8),dt.datetime(2023,6,9)]
stations = None # read all stations
collection = None # read latest collection
lvl = 'l1b'

ds = pyrnet.read_thredds(
    dates=dates,
    campaign=campaign,
    stations=stations,
    collection=collection,
    lvl=lvl
)
```

```{note}
`stations` can also be a list of single station numbers, or just one single station number.
 To select the specific data required.
```

## ... decide which data level to use?
Data processing levels of PyrNet are *l1a*, *l1b*.
Refer to {numref}`tab-datalvl` for an overview.


## ... keep full resolution of geo coordinates (processing raw data)?
From the processing step l1a -> l1b the decision whether to keep lat and lon coordinates time resolution (and interpolate to measurement sample), or average them over the maintenance interval has to be made.

In the configuration file provided for ```pyrnet process l1b -c``` the keyword "*average_latlon*" switches the behaviour.
The default is *True*, so the coordinates are averaged over the maintenance period.