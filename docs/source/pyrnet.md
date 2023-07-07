# Introduction

# Technical Description

## Logger
For data logging, PyrNet uses the [Logomatic v2 Serial SD Datalogger](https://www.sparkfun.com/products/retired/10216) (Sparkfun WIG-10216).


# Setup

# Maintenance

## Quality Flagging

(data-processing)=
# Data Processing

(data-level)=
## Data level
The processed PyrNet data is structured in data levels.
Starting from the raw (logger) data in ASCII fromat, the processed data is stored in netCDF4 files

```{list-table} PyrNet data level
:header-rows: 1
:name: tab-datalvl

* - Level
  - Description
  - Resolution
  - Meta data
* - l1a
  - Maintenance period data; Uncalibrated raw readings from logger
  - Full time resolution (10Hz ADC, 1Hz GPS) 
  - Maintenance quality flags; Campaign meta data
* - l1b
  - Daily data; Calibrated;  Resampled to desired resolution
  - User configured time resolution (default: 1Hz); GPS-synchronized datetime index
  - Same as l1a, added calibration information
```



