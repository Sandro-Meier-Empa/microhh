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

import os
from microhh_py import microhh_tools as mht
import argparse
import collections
import glob
import numpy as np
from multiprocessing import Pool, set_start_method
import platform

if platform.system() == 'Darwin':
    try:
        set_start_method("fork")
    except RuntimeError:
        pass


def convert_to_nc_worker(args_tuple):
    """
    Worker function for multiprocessing.
    Unpacks arguments and calls the actual conversion.
    """
    variables, config = args_tuple

    # Extract config for local scope
    modes = config["modes"]
    starttime = config["starttime"]
    iotimeprec = config["iotimeprec"]
    indexes = config["indexes"]
    itot, jtot, ktot = config["itot"], config["jtot"], config["ktot"]
    grid = config["grid"]
    precision = config["precision"]
    compression = config["compression"]
    overwrite = config["overwrite"]
    endtime = config["endtime"]
    sampletime = config["sampletime"]

    # Loop over the different variables and crosssections
    for variable in variables:
        for mode in modes:
            try:
                otime = int(round(starttime / 10**iotimeprec))

                if os.path.isfile("{0}.xy.000.{1:07d}".format(variable, otime)):
                    if mode != 'xy':
                        continue
                    at_surface = True
                else:
                    at_surface = False

                filename = "{0}.{1}.nc".format(variable, mode)
                if os.path.isfile(filename):
                    if overwrite:
                        os.remove(filename)
                        print("Overwriting %s" % filename)
                    else:
                        print("%s already exists. Skipping..." % filename)
                        continue

                halflevel = '000'
                if not at_surface:
                    if indexes is None:
                        indexes_local, halflevel = mht.get_cross_indices(variable, mode)
                    else:
                        indexes_local = indexes

                        files = glob.glob("{0:}.{1}.*.{2:05d}.{3:07d}".format(
                                variable, mode, indexes_local[0], otime))
                        if len(files) == 0:
                            raise Exception('Cannot find any cross-section')
                        halflevel = files[0].split('.')[-3]

                dim = collections.OrderedDict()
                dim['time'] = []
                dim['z'] = range(ktot)
                dim['y'] = range(jtot)
                dim['x'] = range(itot)

                if at_surface:
                    dim.pop('z')
                    n = itot * jtot
                    indexes_local = [-1]
                elif mode == 'xy':
                    dim.update({'z': []})
                    n = itot * jtot
                elif mode == 'xz':
                    dim.update({'y': []})
                    n = itot * ktot
                elif mode == 'yz':
                    dim.update({'x': []})
                    n = ktot * jtot

                if halflevel[0] == '1':
                    dim['xh'] = dim.pop('x')
                if halflevel[1] == '1':
                    dim['yh'] = dim.pop('y')
                if halflevel[2] == '1':
                    dim['zh'] = dim.pop('z')

                ncfile = mht.Create_ncfile(
                    grid, filename, variable, dim, precision, compression)

                for key, val in dim.items():
                    if key == 'time':
                        continue
                    elif val == []:
                        ncfile.dimvar[key][:] = grid.dim[key][indexes_local]

                for t, time in enumerate(np.arange(starttime, endtime + sampletime, sampletime)):
                    for k in range(len(indexes_local)):
                        index = indexes_local[k]
                        otime = int(round((time) / 10**iotimeprec))

                        if at_surface:
                            f_in = "{0}.{1}.{2}.{3:07d}".format(
                                variable, mode, halflevel, otime)
                        else:
                            f_in = "{0:}.{1}.{2}.{3:05d}.{4:07d}".format(
                                variable, mode, halflevel, index, otime)
                        try:
                            fin = mht.Read_binary(grid, f_in)
                        except Exception as ex:
                            print(ex)
                            break

                        print(
                            "Processing %8s, time=%7i, index=%4i"
                            % (variable, otime, index)
                        )

                        ncfile.dimvar['time'][t] = time

                        if at_surface:
                            ncfile.var[t, :, :] = fin.read(n)
                        elif mode == 'xy':
                            ncfile.var[t, k, :, :] = fin.read(n)
                        elif mode == 'xz':
                            ncfile.var[t, :, k, :] = fin.read(n)
                        elif mode == 'yz':
                            ncfile.var[t, :, :, k] = fin.read(n)

                        fin.close()
                ncfile.close()

            except Exception as ex:
                print(ex)
                print("Failed to create %s" % filename)


def run_conversion(
    filename,
    directory=None,
    modes=None,
    variables=None,
    indexes=None,
    starttime=None,
    endtime=None,
    sampletime=None,
    precision="single",
    nprocs=1,
    compression=True,
    order=None,
    overwrite=False,
):
    """
    Run the MicroHH cross-section binary -> NetCDF conversion.

    Parameters mirror the CLI flags:
    - directory (str): working directory (-d/--directory)
    - filename (str): namelist ini file (-f/--filename)
    - vars (list[str] or str): variable names (-v/--vars)
    - modes (list[str]): cross-section modes ('xy', 'xz', 'yz') (-m/--modes)
    - indexes (list[int]): indices (-x/--index)
    - precision ('single'|'double'|None): NetCDF precision (-p/--precision)
    - order (int|None): spatial order 2 or 4 (-o/--order)
    - starttime, endtime, sampletime (float|None): time controls (-t0/-t1/-tstep)
    - nocompression (bool): disable compression (-c/--nocompression)
    - nprocs (int|None): number of worker processes (-n/--nprocs)
    """
    cross_modes = ["xy", "xz", "yz"]

    # 1) Working directory
    if directory is not None:
        os.chdir(directory)

    # 2) Namelist
    nl = mht.Read_namelist(filename)

    itot = nl["grid"]["itot"]
    jtot = nl["grid"]["jtot"]
    ktot = nl["grid"]["ktot"]

    # 3) Time settings
    starttime = float(starttime) if starttime is not None else nl["time"]["starttime"]
    endtime = float(endtime) if endtime is not None else nl["time"]["endtime"]
    sampletime = (
        float(sampletime) if sampletime is not None else nl["cross"]["sampletime"]
    )

    # 4) Modes
    if modes is None:
        modes = list(nl["cross"].keys() & cross_modes)
        if "xy" not in modes:
            for v in np.atleast_1d(nl["cross"]["crosslist"]):
                if "path" in v:
                    modes.append("xy")
                    break

    # 5) Variables
    iotimeprec = nl["time"].get("iotimeprec", 0.0)

    if variables is None:
        variables = nl["cross"]["crosslist"]

    if not isinstance(variables, list):
        variables = [variables]

    if order is None:
        order = nl["grid"].get("swspatialorder", 2)

    grid = mht.Read_grid(itot, jtot, ktot, order=order)

    # Bundle configuration for workers
    config = {
        "modes": modes,
        "starttime": starttime,
        "endtime": endtime,
        "sampletime": sampletime,
        "iotimeprec": iotimeprec,
        "indexes": indexes,
        "itot": itot,
        "jtot": jtot,
        "ktot": ktot,
        "grid": grid,
        "precision": precision,
        "compression": compression,
        "overwrite": overwrite,
    }

    nprocs = min(nprocs, len(variables))
    chunks = [(variables[i::nprocs], config) for i in range(nprocs)]

    with Pool(processes=nprocs) as pool:
        pool.map(convert_to_nc_worker, chunks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MicroHH binary cross-sections to netCDF4 files."
    )
    parser.add_argument(
        "-m",
        "--modes",
        nargs="*",
        help="mode of the cross section",
        choices=["xy", "xz", "yz"],
    )
    parser.add_argument("-f", "--filename", help="ini file name", required=True)
    parser.add_argument("-d", "--directory", help="directory")
    parser.add_argument("-v", "--vars", nargs="*", help="variable names")
    parser.add_argument("-x", "--index", nargs="*", help="indices", type=int)
    parser.add_argument("-t0", "--starttime", help="first time step to be parsed")
    parser.add_argument("-t1", "--endtime", help="last time step to be parsed")
    parser.add_argument("-tstep", "--sampletime", help="time interval to be parsed")
    parser.add_argument(
        "-p",
        "--precision",
        help="precision",
        choices=["single", "double"],
        default="single",
    )
    parser.add_argument(
        "-n", "--nprocs", help="Number of processes", type=int, default=1
    )
    parser.add_argument(
        "-c",
        "--nocompression",
        help="do not compress the netcdf file",
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        help="overwrite existing output netcdf files",
        action="store_true",
    )
    parser.add_argument("-o", "--order", help="order", choices=[2, 4], type=int)

    args = parser.parse_args()

    run_conversion(
        filename=args.filename,
        directory=args.directory,
        modes=args.modes,
        variables=args.vars,
        indexes=args.index,
        starttime=args.starttime,
        endtime=args.endtime,
        sampletime=args.sampletime,
        precision=args.precision,
        nprocs=args.nprocs,
        compression=not args.nocompression,
        order=args.order,
        overwrite=args.overwrite,
    )
