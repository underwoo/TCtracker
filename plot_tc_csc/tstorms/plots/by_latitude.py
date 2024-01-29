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

"""Tropical Storm Timeseries Plot Generator

This module generates the Tropical Storm Timeseries plot using tracking data
generated from a GCM.
"""
import argparse
import os
import shutil
import tempfile
import matplotlib.pyplot
import matplotlib.ticker

from .. import argparse as tsargparse
from ..ori import ori

__all__ = [
]

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-o",
                           help="Directory where plots will be stored",
                           metavar="outDir",
                           dest="outDir",
                           default=os.getcwd(),
                           type=tsargparse.absPath,
                           action=tsargparse.createDir)
    argparser.add_argument("-H",
                           help="Indicates plot is number of hurricanes",
                           dest="do_hur",
                           action='store_true')
    argparser.add_argument("--fraction",
                           help="Plot of fraction of storms instead of "
                                "per-year average",
                           metafar="frac",
                           dest="frac",
                           action='store_true')
    argparser.add_argument("inDir",
                           help="Directory where tropical storm data are " +
                                "available",
                           metavar="inDir",
                           type=tsargparse.absPath,
                           action=tsargparse.dirExists)
    argparser.add_argument("obsDir",
                           help="Compare to observations in directory",
                           metavar="obsDir",
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
    argparser.add_argument("expName",
                           help="Experiment name used in plots",
                           metavar="expName",
                           type=str)
    args = argparser.parse_args()

    storm_type = 'Tropical Storm'
    if (args.do_hur):
        storm_type = 'Hurricane (CAT. 1-5)'
    year_range = f"{args.beg_year:04d}-{args.end_year:04d}"

    # obs axis
    obs = ori(args.obsDir, args.beg_year, args.end_year, 'obs')
    if args.frac:
        obs_yaxis = obs.storm_count.sum(0)/obs.storm_count.sum()
    else:
        obs_yaxis = obs.storm_count.sum(0)/(obs.end_year-obs.start_year+1)
    obs_xaxis = [x * obs.DLAT + obs.SLAT + 0.5*obs.DLAT
                 for x in range(len(obs_yaxis))]
    # Model axis
    model = ori(args.inDir, args.beg_year, args.end_year, 'model')
    if args.frac:
        model_yaxis = model.storm_count.sum(0)/model.storm_count.sum()
    else:
        model_yaxis = \
            model.storm_count.sum(0)/(model.end_year - model.start_year + 1)
    model_xaxis = [x*obs.DLAT + obs.SLAT + 0.5*obs.DLAT
                   for x in range(len(obs_yaxis))]

    # Configure the plot
    fig, ax = matplotlib.pyplot.subplots()
    fig.suptitle(f"{storm_type} Count By Latitude ({year_range})")
    obs, = ax.plot(obs_xaxis, obs_yaxis, label="obs")
    model, = ax.plot(model_xaxis, model_yaxis, label=args.expName)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    ax.set_xlim(-40, 40)
    ax.set_xlabel(u"latitude (\N{DEGREE SIGN}N)")
    ax.set_ylim(0)
    ax.set_ylabel(f"no. {storm_type} per year")
    ax.legend(handles=[obs, model])
    ax.grid(True)
    ax.axhline(y=0)

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)

        plot_filename = "by_latitude.pdf"
        matplotlib.pyplot.savefig(plot_filename)
        shutil.copyfile(plot_filename,
                        os.path.join(args.outDir, plot_filename))
        print(f"Plot stored in '{os.path.join(args.outDir, plot_filename)}'")
