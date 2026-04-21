#
#  MicroHH
#  Copyright (c) 2011-2024 Chiel van Heerwaarden
#  Copyright (c) 2011-2024 Thijs Heus
#  Copyright (c) 2014-2024 Bart van Stratum
#
#  This file is part of MicroHH
#
#  MicroHH is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MicroHH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
#

import microhh_py.microhh_tools as mht
import argparse
import os
import glob
import struct
import numpy as np
from multiprocessing import Pool, set_start_method
import platform

if platform.system() == "Darwin":
    try:
        set_start_method("fork")
    except RuntimeError:
        pass


def convert_to_nc_worker(args_tuple):
    variables, config = args_tuple

    itot = config["itot"]
    jtot = config["jtot"]
    kmax = config["kmax"]
    starttime = config["starttime"]
    endtime = config["endtime"]
    sampletime = config["sampletime"]
    iotimeprec = config["iotimeprec"]
    doubledump = config["doubledump"]
    grid = config["grid"]
    precision = config["precision"]
    perslice = config["perslice"]
    compression = config["compression"]
    overwrite = config["overwrite"]

    half_level_vars = ["w", "lflx", "sflx"]
    for variable in variables:
        filename = "{0}.nc".format(variable)

        if os.path.isfile(filename):
            if overwrite:
                os.remove(filename)
                print("Overwriting %s" % filename)
            else:
                print("%s already exists. Skipping..." % filename)
                continue

        dim = {"time": [], "z": range(kmax), "y": range(jtot), "x": range(itot)}
        if variable == "u":
            dim["xh"] = dim.pop("x")
        if variable == "v":
            dim["yh"] = dim.pop("y")
        if variable in half_level_vars:
            dim["zh"] = dim.pop("z")
        try:

            def convert(otime, tout):
                f_in = "{0:}.{1:07d}".format(variable, otime)
                try:
                    fin = mht.Read_binary(grid, f_in)
                except Exception as ex:
                    print(ex)
                    return
                # raise Exception(
                #         'Stopping: cannot find file {}'.format(f_in))
                print("Processing %8s, time=%7i" % (variable, otime))
                ncfile.dimvar["time"][tout] = otime * 10**iotimeprec
                if perslice:
                    for k in range(kmax):
                        ncfile.var[tout, k, :, :] = fin.read(itot * jtot)
                else:
                    ncfile.var[tout, :, :, :] = fin.read(itot * jtot * kmax)

                fin.close()

            ncfile = mht.Create_ncfile(
                grid, filename, variable, dim, precision, compression
            )
            # Loop through the files and read 3d field
            tout = 0
            for t, time in enumerate(
                np.arange(starttime, endtime + sampletime, sampletime)
            ):
                otime = int(round(time / 10**iotimeprec))
                if doubledump and t > 0:
                    with open("time.{0:07d}".format(otime), "rb") as file_handle:
                        timedata = struct.unpack("=QQi", file_handle.read())
                    otime2 = int(
                        round(
                            (timedata[0] - timedata[1]) * 10 ** (-iotimeprec - 9) - 0.5
                        )
                    )
                    convert(otime2, tout)
                    tout += 1

                try:
                    convert(otime, tout)
                except Exception as ex:
                    print(ex)
                    break
                tout += 1
            ncfile.close()
        except Exception as ex:
            print(ex)
            print("Failed to create %s" % filename)


def run_conversion(
    filename,
    directory=None,
    variables=None,
    precision="single",
    order=None,
    starttime=None,
    endtime=None,
    sampletime=None,
    perslice=False,
    compression=True,
    kmax=None,
    nprocs=None,
    overwrite=False,
):
    """
    Run the MicroHH 3D binary -> NetCDF conversion.

    Parameters mirror the CLI flags:
    - directory (str): working directory (-d/--directory)
    - filename (str): namelist ini file (-f/--filename)
    - vars (list[str] or str): variable names (-v/--vars)
    - precision ('single'|'double'|None): NetCDF precision (-p/--precision)
    - order (int|None): spatial order 2 or 4 (-o/--order)
    - starttime, endtime, sampletime (float|None): time controls (-t0/-t1/-tstep)
    - perslice (bool): per-slice IO (-s/--perslice)
    - nocompression (bool): disable compression (-c/--nocompression)
    - kmax (int|None): reduce vertical extent (-kmax/--kmax)
    - nprocs (int|None): number of worker processes (-n/--nprocs)
    """
    # 1) Working directory
    if directory is not None:
        os.chdir(directory)

    nl = mht.Read_namelist(filename)

    itot = nl["grid"]["itot"]
    jtot = nl["grid"]["jtot"]
    ktot = nl["grid"]["ktot"]

    kmax = ktot if kmax is None else min(kmax, ktot)
    starttime = float(starttime) if starttime is not None else nl["time"]["starttime"]
    endtime = float(endtime) if endtime is not None else nl["time"]["endtime"]
    sampletime = (
        float(sampletime) if sampletime is not None else nl["dump"]["sampletime"]
    )

    doubledump = nl.get("dump", default={}).get("swdoubledump", 0) == 1
    iotimeprec = nl.get("time", default={}).get("iotimeprec", 0.0)

    if variables is None:
        variables = nl["dump"]["dumplist"]
    if not isinstance(variables, list):
        variables = [variables]
    if len(variables) == 0:
        return

    if order is None:
        order = nl["grid"].get("swspatialorder", 2)

    for time in np.arange(starttime, endtime, sampletime):
        otime = int(round(time / 10**iotimeprec))
        if not glob.glob("*.{0:07d}".format(otime)):
            endtime = time - sampletime
            break

    grid = mht.Read_grid(itot, jtot, ktot, order=order)

    if kmax < ktot:
        grid.dim["z"] = grid.dim["z"][:kmax]
        grid.dim["zh"] = grid.dim["zh"][: kmax + 1]

    config = {
        "itot": itot,
        "jtot": jtot,
        "kmax": kmax,
        "starttime": starttime,
        "endtime": endtime,
        "sampletime": sampletime,
        "iotimeprec": iotimeprec,
        "doubledump": doubledump,
        "grid": grid,
        "precision": precision,
        "perslice": perslice,
        "compression": compression,
        "overwrite": overwrite,
    }

    nprocs = max(nprocs or 1, len(variables))
    chunks = [(variables[i::nprocs], config) for i in range(nprocs)]

    with Pool(processes=nprocs) as pool:
        pool.map(convert_to_nc_worker, chunks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MicroHH 3D binary to netCDF4 files."
    )
    parser.add_argument("-d", "--directory", help="directory")
    parser.add_argument("-f", "--filename", help="ini file name", required=True)
    parser.add_argument("-v", "--vars", nargs="*", help="variable names")
    parser.add_argument(
        "-p",
        "--precision",
        help="precision",
        choices=["single", "double"],
        default="single",
    )
    parser.add_argument("-o", "--order", help="order", choices=[2, 4], type=int)
    parser.add_argument(
        "-t0", "--starttime", help="first time step to be parsed", type=float
    )
    parser.add_argument(
        "-t1", "--endtime", help="last time step to be parsed", type=float
    )
    parser.add_argument(
        "-tstep", "--sampletime", help="time interval to be parsed", type=float
    )
    parser.add_argument(
        "-s", "--perslice", help="read/write per horizontal slice", action="store_true"
    )
    parser.add_argument(
        "-c",
        "--nocompression",
        help="do not compress the netcdf file",
        action="store_true",
    )
    parser.add_argument(
        "-kmax", "--kmax", help="reduce vertical extent 3D files", type=int
    )
    parser.add_argument(
        "-n", "--nprocs", help="Number of processes", type=int, default=None
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        help="overwrite existing output netcdf files",
        action="store_true",
    )

    args = parser.parse_args()

    run_conversion(
        filename=args.filename,
        directory=args.directory,
        variables=args.vars,
        precision=args.precision,
        order=args.order,
        starttime=args.starttime,
        endtime=args.endtime,
        sampletime=args.sampletime,
        perslice=args.perslice,
        compression=not args.nocompression,
        kmax=args.kmax,
        nprocs=args.nprocs,
        overwrite=args.overwrite,
    )
