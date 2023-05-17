import re
import os.path

import click
import numpy as np
import pandas as pd
import pkg_resources as pkg_res

from . import utils
from . import logger
from . import reports

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
              multiple=True,
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

    # read default config
    fn = pkg_res.resource_filename(__package__, "share/pyrnet_config.json")
    cfg = utils.read_json(fn)

    # update config with input config files
    for fn in config:
        cfg.update(utils.read_json(fn))

    # add default cfmeta if user doesn't specify
    if not cfg['cfjson']:
        fn_cfjson = pkg_res.resource_filename("pyrnet", "share/pyrnet_cfmeta_l1b.json")
        cfg.update({"cfjson": fn_cfjson})

    # filename parser
    parse = re.compile(cfg['filename_parser'])

    # parse maintenance reports
    if report=="" or report=="online":
        df_report = reports.get_responses(fn=None, online=cfg["online"])
    else:
        df_report = reports.get_responses(fn=report)
    report = reports.parse_report(df_report,
                                  date_of_maintenance=np.datetime64(date_of_maintenance))

    with click.progressbar(input_files,label='Processing') as files:
        for fn in files:
            filepath = os.path.abspath(fn)
            filename = os.path.basename(filepath)

            m = parse.match(filename)
            try:
                stationid = int(m.group('ID'))
            except:
                raise ValueError(f"Could not find station id in filename {filename} using regex {config['filename_parser']}.")

            ds = logger.read_logger(
                fname=fn,
                station=stationid,
                bins=cfg['bins'],
                date_of_measure=np.datetime64(cfg['date_of_measure']),
                report=report,
                config=cfg,
                global_attrs=cfg['global_attrs']
            )

            outfile = os.path.join(output_path, cfg['output'])
            outfile = outfile.format_map(
                dict(
                    startdt=pd.to_datetime(ds.time.values[0]),
                    enddt=pd.to_datetime(ds.time.values[-1]),
                    campaign=cfg['campaign'],
                    station=stationid,
                    level='l1b',
                    collection=cfg['collection'],
                    sfx="nc"
                )
            )
            ds.to_netcdf(outfile)


@click.command("l1b")
def process_l1b():
    print("Process l1b")


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