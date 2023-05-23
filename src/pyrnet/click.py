import re
import os.path

import click
import numpy as np
import pandas as pd
import xarray as xr
import pkg_resources as pkg_res

from . import pyrnet
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
            # print(f"{filename}",flush=True)
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

            outfile = os.path.join(output_path, cfg['output_l1a'])
            outfile = outfile.format_map(
                dict(
                    startdt=pd.to_datetime(ds.time.values[0]),
                    enddt=pd.to_datetime(ds.time.values[-1]),
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
              multiple=True,
              help="Specify config files with override the default config.")
@click.option("--calibration",
              nargs=1,
              help="Path to the calibration file. If not specified, pyrnet_calibration.json from pyrnet/share/ is used.")
@click.option("--mapping",
              nargs=1,
              help="Path to the calibration file. If not specified, pyrnet_station_map.json from pyrnet/share/ is used.")
@click.option("--radflux_varname",
              nargs=2,
              default=["ghi", "gti"],
              help="Dataset variable name of radiation flux. The default is ['ghi','gti'].")
def process_l1b(input_files: list[str],
                output_path: str,
                config:list[str],
                calibration: str|None,
                mapping:str|None,
                radflux_varname:list[str]):

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

    with click.progressbar(input_files,label='Processing') as files:
        for fn in files:
            filepath = os.path.abspath(fn)
            filename = os.path.basename(filepath)
            # print(f"{filename}", flush=True)
            ds = xr.open_dataset(filepath)
            # check correct file
            if ds.processing_level != "l1a":
                raise ValueError(f"{fn} is not a l1a file.")

            box = ds.station.values[0]
            boxnumber, serial, cfac = pyrnet.meta_lookup(
                ds.time.values[0],
                box=box,
                cfile=calibration,
                mapfile=mapping,
            )
            # print(f"Box:{boxnumber}, Serials: {serial}, Calibration: {cfac}",flush=True)

            # calibrate radiation flux with gain=300
            for i,radflx in enumerate(radflux_varname):
                ds[radflx].values = ds[radflx].values*1e6/(cfac[i]* 300) # V -> Wm-2
                ds[radflx].attrs['units'] = "W m-2",
                ds[radflx].encoding.update({
                    'scale_factor': ds[radflx].encoding['scale_factor']*1e6/(cfac[i]* 300)
                })


            ds.attrs["processing_level"] = 'l1b'
            now = pd.to_datetime(np.datetime64("now"))
            ds.attrs["history"] = ds.history + f"{now.isoformat()}: Generated level l1a  by pyrnet version XX; " #TODO add version

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
                dsd.to_netcdf(outfile)

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