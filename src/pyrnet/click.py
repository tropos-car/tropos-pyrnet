import re
import os.path

import click
import numpy as np
import pandas as pd
import xarray as xr
import logging
from collections.abc import Iterable
from toolz import merge_with, assoc_in
import parse

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
    if report is None:
        df_report = None
    elif report=="online":
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

            outfile = os.path.join(
                output_path,
                pyrdata.get_fname(ds, freq="10Hz", timevar="gpstime", sfx="nc", config=cfg)
            )

            pyrdata.to_netcdf(ds, outfile, timevar="gpstime")
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
            if ds is None:
                logger.debug(f"{filename} is skipped.")
                continue

            box = int(ds.station.values[0])
            udays = np.unique(ds.time.values.astype("datetime64[D]"))
            for day in udays:
                day = pd.to_datetime(day)
                logging.info(f"process day {day:%Y-%m-%d}")
                dsd = ds.sel(time=f"{day:%Y-%m-%d}")

                outfile = os.path.join(
                    output_path,
                    pyrdata.get_fname(dsd, period="P1D", freq=cfg["l1bfreq"], timevar="time", sfx="nc", config=cfg)
                )

                pyrdata.to_netcdf_l1b(dsd, fname=outfile, freq=cfg["l1bfreq"])
                logging.info(f"l1b saved to {outfile}")

@click.command("l1b_network")
@click.argument("input_files", nargs=-1)
@click.argument("output_path", nargs=1)
@click.option("--config","-c",
              nargs=1,
              help="Specify config files with override the default config.")
def process_l1b_network(input_files: list[str],
                output_path: str,
                config:str):

    if config is not None:
        config = pyrutils.read_json(config)
    cfg = pyrdata.get_config(config)

    # get unique station numbers
    stations = []
    for fn in input_files:
        filename = os.path.basename(fn)
        result = parse.parse(cfg["output"], filename).named
        stations.append(result["station"])
    Nstations = len(np.unique(stations))

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
            if ds is None:
                logger.debug(f"{filename} is skipped.")
                continue

            udays = np.unique(ds.time.values.astype("datetime64[D]"))
            for day in udays:
                day = pd.to_datetime(day)
                logging.info(f"process day {day:%Y-%m-%d}")
                dsd = ds.sel(time=f"{day:%Y-%m-%d}")

                outfile = os.path.join(
                    output_path,
                    pyrdata.get_fname(
                        dsd,
                        period="P1D",
                        kind='n',
                        station=0, #Nstations,
                        freq=cfg["l1bfreq"],
                        timevar="time",
                        sfx="nc",
                        config=cfg
                    )
                )


                pyrdata.to_netcdf_l1b(dsd, fname=outfile, freq=cfg["l1bfreq"])
                logging.info(f"l1b_network saved to {outfile}")


cli.add_command(process)
process.add_command(process_l1a)
process.add_command(process_l1b)
process.add_command(process_l1b_network)

@click.command("merge")
@click.argument("input_files", nargs=-1)
@click.argument("output_file", nargs=1)
@click.option("-f","--freq",nargs=1,help="Sampling frequency for regular time grid. The default is 1s.")
@click.option("-t","--timevar", nargs=1, help="Name of the variable storing the time index. The default is 'time'.")
def merge(input_files, output_file,freq=None,timevar=None):
    def _read_radflux_attrs(ds):
        def _ensure_list(a):
            if (not isinstance(a, Iterable)) or isinstance(a, str):
                return [a]
            else:
                return list(a)

        vattrs = {}
        for var in ['ghi', 'gti']:
            vattrs.update({
                var: {
                    "serial": _ensure_list(ds[var].serial),
                    "calibration_factor": _ensure_list(ds[var].calibration_factor)
                }
            })
            if var == "gti":
                vattrs = assoc_in(vattrs, ["gti", "hangle"],
                                  _ensure_list(ds[var].hangle))
                vattrs = assoc_in(vattrs, ["gti", "vangle"],
                                  _ensure_list(ds[var].vangle))

        return vattrs

    if timevar is None:
        timevar = "time"
    if freq is None:
        freq = '1s'

    with click.progressbar(input_files, label='Merging') as files:

        for i, fn in enumerate(files):
            dst = xr.open_dataset(fn)
            # unify time dimension to speed up merging
            date = dst.time.values[0].astype("datetime64[D]")
            timeidx = pd.date_range(date, date + np.timedelta64(1, 'D'), freq=freq, inclusive='left')
            dst = dst.reindex(time=timeidx, method='nearest', tolerance=np.timedelta64(1, 'ms'))

            # add gti for single stations
            if "gti" not in dst:
                dst = dst.assign({
                    "gti": (dst.ghi.dims, np.full(dst.ghi.values.shape, np.nan)),
                    "qc_flag_gti": (dst.qc_flag_ghi.dims, np.full(dst.qc_flag_ghi.values.shape, 0)),
                    "maintenance_flag_gti": (
                    dst.maintenance_flag_ghi.dims, np.full(dst.maintenance_flag_ghi.values.shape, 0))
                })
                dst.gti.attrs.update({
                    "serial":"",
                    "calibration_factor": 0,
                    "vangle": 0,
                    "hangle": 0
                })

            if i==0:
                ds = dst.copy()
                vattrs_radflx = _read_radflux_attrs(ds)
                continue

            st = dst.station.values[0]
            if st not in ds.station.values:
                vattrs_temp = _read_radflux_attrs(dst)
                vattrs_radflx.update({
                    "ghi": merge_with(lambda x: [*x[0],*x[1]], (vattrs_radflx['ghi'], vattrs_temp['ghi'])),
                    "gti": merge_with(lambda x: [*x[0],*x[1]], (vattrs_radflx['gti'], vattrs_temp['gti']))
                })
                ds = xr.concat((ds, dst), dim='station')
            else:
                overwrite_vars = [v for v in ds if timevar not in ds[v].dims]
                ds = ds.merge(dst,
                              compat='no_conflicts',
                              overwrite_vars=overwrite_vars)
            dst.close()

    # special treatment for flux variables
    for k in ['ghi', 'gti']:
        if k not in ds:
            continue
        # add concatenated attrs
        ds[k].attrs.update(vattrs_radflx[k])

    ds = pyrdata.add_encoding(ds)

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