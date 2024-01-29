# **********************************************************************
# TCtracker - Tropical Storm Detection
# Copyright (C) 2021, 2023 Frederic Vitart, Joe Sirutis, Ming Zhao,
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

"""
This module contains classes and other helper routines to read
ori_<year> files that contain tropical storm data.

The ori_<year> data has the format:
longitude(float), latitude(float), year(int), month(int), day(int), hour(int)
"""

import tempfile
import os
import shutil
import re
import numpy

from ..config import pkglibexecdir  # pyright: ignore [reportMissingImports]
from .stat_ori_mask import stat_ori as _stat_ori  # pyright: ignore [reportMissingImports] # noqa: E501
from .StormBox import StormBox
from .freq_ori import freq_ori as _freq_ori  # pyright: ignore [reportMissingImports] # noqa: E501

__all__ = ["ori",
           "StormBox"]


class ori():
    """Class to hold data for a group of `ori_<year>` files

    Keyword Arguments:

        - ori_dir -- Directory that contains the `ori_<year>` files.

        - beg_year -- First year of `ori_<year>` data to process.

        - end_year -- Last year of `ori_<year>` data to process.

        - ori_type -- Type of data in directory.  Is typically 'obs' for
            the observational `ori_<year>` files, or 'model' for model
            generated `ori_<year>` files.  Default 'model'

    Attributes:

        - directory -- Holds `ori_dir`

        - start_year -- Holds `beg_year`

        - end_year -- Holds `end_year`

        - type -- Holds `ori_type`

        - storm_count -- Holds the total number of storms per 5x4 lat/lon box

        - stat_file -- Name of the statistics file using `ori_<start_year>`-
            `ori_<end_year>` files.

        - stats -- Dictionary of StormBox with region IDs as key.  Data
            obtained form data in `stat_file`.
    """

    DLON = 5.0  # Longitude region delta
    DLAT = 4.0  # Latitude region delta
    SLAT = -86.0  # Lowest southern latitude

    def __init__(self,
                 ori_dir: str,
                 beg_year: int,
                 end_year: int,
                 ori_type='model'):
        self.directory = ori_dir
        self.start_year = beg_year
        self.end_year = end_year
        self.type = ori_type
        self.stat_file = self._gen_stats()
        self.stats = self._read_stats()

        # Calculated Parameters
        ix = int(360.0/self.DLON)
        jx = int((90-self.SLAT)/self.DLAT)
        # Closest indexes to 40s and 40n
        # j40s = int((-40.0-self.SLAT)/self.DLAT) + 1
        # j40n = jx - j40s + 1

        # Hold the storm frequency map
        self.storm_count = numpy.zeros((ix, jx), dtype=int)

        for year in range(beg_year, end_year + 1):
            # I'm assuming all ori_<year> files are available (even if empty)
            # Will need to do some checking in the future
            fname = os.path.join(self.directory, "ori_{:04d}".format(year))
            with open(fname, 'r') as infile:
                for line in infile:
                    xcyc, ycyc, year, month, day, hour = line.split()
                    # Until I find something better
                    # year, month, day, hour are not used ATT
                    xcyc = float(xcyc)
                    ycyc = float(ycyc)
                    year = int(year)
                    month = int(month)
                    day = int(day)
                    hour = int(hour)

                    # These are indices for the lon/lat frequency arrays
                    # The added 0.5 is to center the values around the
                    # DLAT/DLON divisions
                    i = int(xcyc/self.DLON + 0.5)
                    j = int((ycyc - self.SLAT)/self.DLAT + 0.5)

                    # Deal with the wrap around indices
                    if i < 0:
                        i = ix + i
                    elif i >= ix:
                        i = i - ix
                    self.storm_count[i][j] += 1

    def cat_ori_files(self, fname="ori"):
        """
        Concatinate ori_[YYYY] files into a single `ori` file.  The single
        `ori` file is required for certain plots.
        """

        # Concatenate all `ori_YYYY` files into a single `ori` file
        with open(fname, 'w') as outfile:
            for year in range(self.start_year, self.end_year + 1):
                fname = os.path.join(self.directory, "ori_{:04d}".format(year))
                with open(fname) as infile:
                    outfile.write(infile.read())
        return os.path.realpath(fname)

    def freq_ori(self,
                 do_40ns=True,
                 do_map=True,
                 do_lon=False,
                 do_lat=False,
                 do_latf=False,
                 do_fot=False,
                 traj_in=False):
        """Run the freq_ori Fortran function on the ori data

        This routine concatinates the ori_year files into a single `ori` file,
        and then runs the freq_ori Fortran subroutine on the ori data.  This
        routine will generate a series of files that can be used to generate
        2D plots using Grace.

        Keyword Arguments:

            - do_40ns -- Bound the search between the 40°S and 40°N latitude.
                Only affects `do_map`, `do_lat` and `do_latf`.

            - do_map -- Produce frequency data for full globe.  Output placed
                in file `fmap`.

            - do_lon -- Produce frequency data based on longitude.  Output
                placed in three region files: `flon_gl` (global), `flon_nh`
                (northern hemisphere) and `flon_sh` (southern hemisphere) with
                output:

                <longitude> <frequency>

            - do_lat -- Produce frequency data based on latitude.  Output
                placed in `flat` file with output:

                <latitude> <frequency>

            - do_latf -- Same as `do_lat`, but output is switched:
                <frequency> <latitude>

            - do_fot -- Write frequency as a fraction of total number of
                global storms.  Otherwise, frequency is fraction of storms
                per year.

            - traj_in -- Read in `traj` formated data
        """

        # Ensure the _correct_ ori file is in place
        self.cat_ori_files()
        _freq_ori(do_40ns, do_map, do_lon, do_lat, do_latf, do_fot, traj_in)

    def _gen_stats(self):
        """Generate storm statics from `ori_<year>` files

        Runs the stat_ori Fortran subroutine on the `ori` files, and geneartes
        monthly and annual statistics for different global regions.
        """

        stats_filename = 'stats_{0}_{1:04d}-{2:04d}'.format(self.type,
                                                            self.start_year,
                                                            self.end_year)
        # Remember where we are
        prev_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            # Concatinate all ori_YYYY files into a single ori file
            self.cat_ori_files()
            # Run the ori_stat executable
            _stat_ori(os.path.join(pkglibexecdir, 'imask_2'), False, False)

            shutil.copyfile('stat_mo', os.path.join(prev_cwd, stats_filename))

        os.chdir(prev_cwd)
        return os.path.realpath(stats_filename)

    def _read_stats(self):
        """Read in stat file data

        Open a stat file, read in the contents, and return a dict of
        StormBox's with the Box ID as keys.
        """

        # Regex patterns for box, storms (by year) and stats over the time
        # period for the box.
        __boxhdr = re.compile(r"^ {2}\*{3} Box = {1,2}([A-Z]{1,2})$")
        __yearln = re.compile(r"^ *([0-9]+)((?: *[0-9]+){12}) *([0-9]+)$")
        __statln = re.compile(r"^ *(sum|sprd|mean|std)((?: *[0-9.]+){13})$")

        return_boxes = {}
        with open(self.stat_file, 'r') as f:
            for line in f:
                box_match = __boxhdr.match(line)
                if box_match:
                    box = StormBox(box_match.group(1))
                    next(f, None)  # The next line is the header line.
                    for line2 in f:
                        year_match = __yearln.match(line2)
                        stat_match = __statln.match(line2)
                        if year_match:
                            box.add_storms(year_match.group(1),
                                           year_match.group(2),
                                           year_match.group(3))
                        elif stat_match:
                            box.add_stats(stat_match.group(1),
                                          stat_match.group(2))
                        else:
                            return_boxes[box.id] = box
                            break
        return return_boxes
