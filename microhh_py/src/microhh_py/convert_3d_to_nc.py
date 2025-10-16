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
import time as tm
import numpy as np
from multiprocessing import Pool


def convert_to_nc(variables):
    half_level_vars = ["w", "lflx", "sflx"]

    for variable in variables:
        filename = "{0}.nc".format(variable)
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
                otime = round(time / 10**iotimeprec)
                if doubledump and t > 0:
                    timedata = struct.unpack(
                        "=QQi", open("time.{0:07d}".format(otime), "rb").read()
                    )
                    otime2 = round(
                        (timedata[0] - timedata[1]) * 10 ** (-iotimeprec - 9) - 0.5
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


def run(
    directory=None,
    filename=None,
    vars=None,
    precision=None,
    order=None,
    starttime=None,
    endtime=None,
    sampletime=None,
    perslice=False,
    nocompression=False,
    kmax=None,
    nprocs=None,
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

    # 2) Namelist
    if not filename:
        raise ValueError("filename (namelist ini) must be provided")
    nl = mht.Read_namelist(filename)
    itot = nl["grid"]["itot"]
    jtot = nl["grid"]["jtot"]
    ktot = nl["grid"]["ktot"]
    kmax_local = min(kmax if kmax is not None else ktot, ktot)

    # 3) Time & dump settings
    starttime = starttime if starttime is not None else nl["time"]["starttime"]
    endtime = endtime if endtime is not None else nl["time"]["endtime"]
    sampletime = sampletime if sampletime is not None else nl["dump"]["sampletime"]
    try:
        doubledump = nl["dump"]["swdoubledump"] == 1
    except Exception:
        doubledump = False

    try:
        iotimeprec = nl["time"]["iotimeprec"]
    except KeyError:
        iotimeprec = 0.0

    variables = vars if vars is not None else nl["dump"]["dumplist"]
    if isinstance(variables, str):
        variables = [variables]

    # promote to globals used by convert_to_nc
    globals().update(
        {
            "itot": itot,
            "jtot": jtot,
            "ktot": ktot,
            "kmax": kmax_local,
            "starttime": starttime,
            "endtime": endtime,
            "sampletime": sampletime,
            "doubledump": doubledump,
            "iotimeprec": iotimeprec,
            "precision": precision,
            "perslice": perslice,
            "compression": not nocompression,
        }
    )

    try:
        order = order if order is not None else nl["grid"]["swspatialorder"]
    except KeyError:
        order = 2

    # 4) Truncate endtime to last available dump
    for time in np.arange(starttime, endtime, sampletime):
        otime = int(round(time / 10**iotimeprec))
        if not glob.glob("*.{0:07d}".format(otime)):
            endtime = time - sampletime
            break
    globals().update({"endtime": endtime})

    # 5) Grid
    grid = mht.Read_grid(itot, jtot, ktot, order=order)
    if kmax_local < ktot:
        grid.dim["z"] = grid.dim["z"][:kmax_local]
        grid.dim["zh"] = grid.dim["zh"][: kmax_local + 1]
    globals().update({"grid": grid})

    # 6) Parallel chunks
    nprocs = nprocs if nprocs is not None else len(variables)
    chunks = [variables[i::nprocs] for i in range(max(1, nprocs))]

    # 7) Run
    with Pool(processes=nprocs) as pool:
        for _ in pool.imap_unordered(convert_to_nc, chunks):
            pass  # progress is printed inside convert_to_nc


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description="Convert MicroHH 3D binary to netCDF4 files."
    )
    p.add_argument("-d", "--directory", help="directory")
    p.add_argument("-f", "--filename", help="ini file name")
    p.add_argument("-v", "--vars", nargs="*", help="variable names")
    p.add_argument("-p", "--precision", choices=["single", "double"])
    p.add_argument("-o", "--order", choices=[2, 4], type=int)
    p.add_argument(
        "-t0", "--starttime", type=float, help="first time step to be parsed"
    )
    p.add_argument("-t1", "--endtime", type=float, help="last time step to be parsed")
    p.add_argument(
        "-tstep", "--sampletime", type=float, help="time interval to be parsed"
    )
    p.add_argument(
        "-s", "--perslice", action="store_true", help="read/write per horizontal slice"
    )
    p.add_argument(
        "-c",
        "--nocompression",
        action="store_true",
        help="do not compress the netcdf file",
    )
    p.add_argument("-kmax", "--kmax", type=int, help="reduce vertical extent 3D files")
    p.add_argument("-n", "--nprocs", type=int, help="Number of processes")
    return p


def main():
    args = _build_arg_parser().parse_args()
    run(
        directory=args.directory,
        filename=args.filename,
        vars=args.vars,
        precision=args.precision,
        order=args.order,
        starttime=args.starttime,
        endtime=args.endtime,
        sampletime=args.sampletime,
        perslice=args.perslice,
        nocompression=args.nocompression,
        kmax=args.kmax,
        nprocs=args.nprocs,
    )


if __name__ == "__main__":
    main()
