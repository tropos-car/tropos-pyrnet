import re
import os.path

import click
import numpy as np
import pandas as pd
import xarray as xr
import logging
from toolz import merge_with

from . import pyrnet
from . import data as pyrdata
from . import utils as pyrutils
from . import reports as pyrreports


# logging setup
logging.basicConfig(
    filename='pyrnet.log',
    encoding='utf-8',
    level=logging.DEBUG,
    format='%(asctime)s %(name)s %(levelname)s:%(message)s'
)
logger = logging.getLogger(__name__)

@click.group("pyrnet")
def cli():
    pass

@click.group("process")
def process():
    print("Process")

@click.command("l1a")
@click.argument("input_files", nargs=-1)
@click.argument("output_path", nargs=1)
@click.option("--config","-c",
              nargs=1,
              help="Specify config files with override the default config.")
@click.option("--report","-r",
              default="",
              help="Specify the maintenance report file. If empty or 'online' it attempts to request it online.")
@click.option("--date_of_maintenance",
              help="Specify date of maintenance as datetime64 string ('YYYY-MM-DD'). If not specified, try to retrieve from data.")
def process_l1a(input_files,
                output_path,
                config,
                report,
                date_of_maintenance):
    if config is not None:
        config = pyrutils.read_json(config)
    cfg = pyrdata.get_config(config)

    # filename parser
    parse = re.compile(cfg['filename_parser'])

    # parse maintenance reports
    if report=="" or report=="online":
        df_report = pyrreports.get_responses(fn=None, online=cfg["online"])
    else:
        df_report = pyrreports.get_responses(fn=report)

    if date_of_maintenance is None:
        report = df_report
    else:
        report = pyrreports.parse_report(df_report,
                                  date_of_maintenance=np.datetime64(date_of_maintenance))

    with click.progressbar(input_files,label='Processing') as files:
        for fn in files:
            filepath = os.path.abspath(fn)
            filename = os.path.basename(filepath)
            logging.info(f"start raw->l1a: {filename}")

            m = parse.match(filename)
            try:
                stationid = int(m.group('ID'))
            except:
                raise ValueError(f"Could not find station id in filename {filename} using regex {config['filename_parser']}.")
            logging.info(f"found station number {stationid}")

            ds = pyrdata.to_l1a(
                fname=fn,
                station=stationid,
                date_of_measure=np.datetime64(cfg['date_of_measure']),
                report=report,
                config=cfg,
                global_attrs=cfg['global_attrs']
            )
            if ds is None:
                logging.warning(f"Skip {filename}.")
                continue

            outfile = os.path.join(output_path, cfg['output_l1a'])
            outfile = outfile.format_map(
                dict(
                    startdt=pd.to_datetime(ds.gpstime.values[0]),
                    enddt=pd.to_datetime(ds.gpstime.values[-1]),
                    campaign=cfg['campaign'],
                    station=stationid,
                    collection=cfg['collection'],
                    sfx="nc"
                )
            )
            ds.to_netcdf(outfile, encoding={'gpstime':{'dtype':'float64'}})
            logging.info(f"l1a saved to {outfile}")


@click.command("l1b")
@click.argument("input_files", nargs=-1)
@click.argument("output_path", nargs=1)
@click.option("--config","-c",
              nargs=1,
              help="Specify config files with override the default config.")
def process_l1b(input_files: list[str],
                output_path: str,
                config:str):

    if config is not None:
        config = pyrutils.read_json(config)
    cfg = pyrdata.get_config(config)

    with click.progressbar(input_files,label='Processing') as files:
        for fn in files:
            filepath = os.path.abspath(fn)
            filename = os.path.basename(filepath)
            logging.info(f"start l1a->l1b: {filename}")

            ds = pyrdata.to_l1b(
                filepath,
                config=config,
                global_attrs=cfg['global_attrs']
            )
            box = ds.station.values[0]
            udays = np.unique(ds.time.values.astype("datetime64[D]"))
            for day in udays:
                day = pd.to_datetime(day)
                logging.info(f"process day {day:%Y-%m-%d}")
                dsd = ds.sel(time=f"{day:%Y-%m-%d}")
                outfile = os.path.join(output_path, cfg['output_l1b'])
                outfile = outfile.format_map(
                    dict(
                        dt=day,
                        campaign=cfg['campaign'],
                        station=box,
                        collection=cfg['collection'],
                        sfx="nc"
                    )
                )
                pyrdata.to_netcdf(dsd,outfile)
                logging.info(f"l1b saved to {outfile}")

cli.add_command(process)
process.add_command(process_l1a)
process.add_command(process_l1b)

@click.command("merge")
@click.argument("input_files", nargs=-1)
@click.argument("output_file", nargs=1)
def merge(input_files, output_file):
    def _read_radflux_attrs(ds):
        vattrs = {}
        for var in ['ghi','gti']:
            if var not in ds:
                vattrs.update({
                    var: {
                        "serial": [""],
                        "calibration_factor": [np.nan]
                    }
                })
                continue
            vattrs.update({
                var: {
                    "serial": [ds[var].serial],
                    "calibration_factor": [ds[var].calibration_factor]
                }
            })
        if "gti" in ds:
            vattrs.update({
                var: {
                    "hangle": [ds[var].hangle],
                    "vangle": [ds[var].vangle]
                }
            })
        return vattrs

    with click.progressbar(input_files, label='Merging') as files:

        for i, fn in enumerate(files):
            if i==0:
                ds = xr.open_dataset(fn)
                vattrs = _read_radflux_attrs(ds)
                continue

            dst = xr.open_dataset(fn)
            st = dst.station.values[0]
            if st not in ds.station.values:
                vattrs_temp = _read_radflux_attrs(dst)
                vattrs.update({
                    "ghi": merge_with(lambda x: [*x[0],*x[1]],(vattrs['ghi'],vattrs_temp['ghi'])),
                    "gti": merge_with(lambda x: [*x[0], *x[1]], (vattrs['gti'], vattrs_temp['gti']))
                })
            ds = xr.merge((ds, dst))

    # prepare netcdf encoding
    gattrs, vattrs, vencode = pyrdata.get_cfmeta()

    # add encoding to Dataset
    for k, v in vencode.items():
        if k not in ds.keys():
            continue
        ds[k].encoding = v

    # special treatment for flux variables
    for k in ['ghi', 'gti']:
        if k not in ds:
            continue
        # add concatenated attrs
        ds[k].attrs.update(vattrs[k])
        # add encoding
        dtype = ds[k].encoding['dtype']
        int_limit = np.iinfo(dtype).max
        valid_range = [0, int_limit - 1]
        scale_factor = 1500. / float(int_limit - 1)
        ds[k].encoding.update({
            "scale_factor": scale_factor,
            "_FillValue": int_limit,
        })
        ds[k].attrs.update({
            "valid_range": valid_range
        })

    ds = pyrdata.stretch_resolution(ds)

    ds["time"].encoding.update({
        "dtype": 'f8',
        "units": f"seconds since {np.datetime_as_string(ds.time.data[0], unit='D')}T00:00Z",
    })
    ds.to_netcdf(output_file)

cli.add_command(merge)



@click.group("convert")
def convert():
    print("Convert")

@click.command("mesor")
def nc2mesor():
    print("convert to mesor.")

cli.add_command(convert)
convert.add_command(nc2mesor)