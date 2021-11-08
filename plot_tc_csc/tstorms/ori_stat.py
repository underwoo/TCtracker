# **********************************************************************
# TCtracker - Tropical Storm Detection
# Copyright (C) 2021 Frederic Vitart, Joe Sirutis, Ming Zhao,
# Kyle Olivo, Keren Rosado and Seth Underwood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
# **********************************************************************

import argparse
import os
import shutil
import tempfile
import subprocess

from . import argparse as tsargparse
from .config import exeext, pkglibexecdir, gracebat


def cat_ori_files(inDir: str, beg_year: int, end_year: int):
    """
    Concatinate ori_[YYYY] files into a single `ori` file.  The single `ori`
    file is required for certain plots.
    """

    # Concatenate all `ori_YYYY` files into a single `ori` file
    with open('ori', 'w') as outfile:
        for year in range(beg_year, end_year + 1):
            fname = os.path.join(inDir, "ori_{:04d}".format(year))
            with open(fname) as infile:
                outfile.write(infile.read())


def run_stats(inDir: str, beg_year: int, end_year: int):
    """
    Generate the tropical storms stat file
    """

    # Concatinate all ori_YYYY files into a single ori file
    cat_ori_files(inDir, beg_year, end_year)
    # Run the ori_stat executable
    subprocess.run([os.path.join(pkglibexecdir, 'stat_ori_mask' + exeext)], input=b'&input /\n')


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-o",
                           metavar="outDir",
                           help="Directory to place output file",
                           dest="outDir",
                           default=os.getcwd(),
                           type=tsargparse.absPath,
                           action=tsargparse.createDir)
    argparser.add_argument("inDir",
                           help="Directory where tropical storm data are available",
                           metavar="inDir",
                           type=tsargparse.absPath,
                           action=tsargparse.dirExists)
    argparser.add_argument("beg_year",
                           help="First year to process",
                           metavar="beg_year",
                           type=int)
    argparser.add_argument("end_year",
                           help="Last year to process",
                           metavar="end_year",
                           type=int)
    argparser.add_argument("statName",
                           help="Name of statistics type (e.g. model, obs)",
                           metavar="name",
                           default="model",
                           nargs='?',
                           type=str)
    args = argparser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        run_stats(args.inDir,
                  args.beg_year,
                  args.end_year)
        shutil.copyfile('stat_mo',
                        os.path.join(args.outDir,
                                     'stats_{0}_{1:04d}-{2:04d}'.format(args.statName,
                                                                        args.beg_year,
                                                                        args.end_year)))