import re
import os.path

import click
import numpy as np
import pandas as pd
import xarray as xr
import pkg_resources as pkg_res

from . import pyrnet
from . import data as pyrdata
from . import utils as pyrutils
from . import reports as pyrreports

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
              default="now",
              help="Specify date of maintenance as datetime64 string ('YYYY-MM-DD'). The default is 'now' - the most recent reports.")
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
    report = pyrreports.parse_report(df_report,
                                  date_of_maintenance=np.datetime64(date_of_maintenance))

    with click.progressbar(input_files,label='Processing') as files:
        for fn in files:
            filepath = os.path.abspath(fn)
            filename = os.path.basename(filepath)
            # print(f"{filename}",flush=True)
            m = parse.match(filename)
            try:
                stationid = int(m.group('ID'))
            except:
                raise ValueError(f"Could not find station id in filename {filename} using regex {config['filename_parser']}.")

            ds = pyrdata.to_l1a(
                fname=fn,
                station=stationid,
                date_of_measure=np.datetime64(cfg['date_of_measure']),
                report=report,
                config=cfg,
                global_attrs=cfg['global_attrs']
            )

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
            ds.to_netcdf(outfile)


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

            ds = pyrdata.to_l1b(
                filepath,
                config=config,
                global_attrs=cfg['global_attrs']
            )
            box = ds.station.values[0]
            udays = np.unique(ds.time.values.astype("datetime64[D]"))
            for day in udays:
                day = pd.to_datetime(day)
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

cli.add_command(process)
process.add_command(process_l1a)
process.add_command(process_l1b)


@click.group("convert")
def convert():
    print("Convert")

@click.command("mesor")
def nc2mesor():
    print("convert to mesor.")

cli.add_command(convert)
convert.add_command(nc2mesor)