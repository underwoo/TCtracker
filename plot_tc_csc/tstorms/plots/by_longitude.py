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

"""Tropical Storm Timeseries Plot Generator

This module generates the Tropical Storm Timeseries plot using tracking data
generated from a GCM.
"""
import argparse
import os
import shutil
import tempfile
import numpy
import matplotlib.pyplot
import matplotlib.ticker

from .. import argparse as tsargparse
from ..ori import ori

__all__ = []


def generate_axes(ori, frac=False):
    """Extract the axis data for lon, global, Northern and Southern
    Hemispheres

    Parameters
    ----------
    ori : Class ori
        The ori class to extract the axis data from
    frac : bool, optional
        Flag to determine if calcuating the y-axis data using the number of
        total storms (True), or the number of years (Default)

    Returns
    -------
    list * 4
        Returns the xaxis, global, N_hemisphere, S_hemisphere
    """
    xaxis = [x*ori.DLON + 0.5*obs.DLON
             for x in range(len(ori.storm_count.shape[0]))]
    nh_mask = numpy.zeros(obs.storm_count.shape, dtype='bool')
    nh_mask[:, obs_xaxis.index(0):] = True
    if frac:
        glo_axis = ori.storm_count.sum(1)/ori.storm_count.sum()
        nhemi_axis = \
            ori.storm_count.sum(1, where=nh_mask)/ori.storm_count.sum()
    else:
        num_years = ori.end_year - ori.start_year + 1
        glo_axis = obs.storm_count.sum(1)/num_years
        nhemi_axis = ori.storm_count.sum(1, where=nh_mask)/num_years
    shemi_axis = glo_axis - nhemi_axis

    return xaxis, glo_axis, nhemi_axis, shemi_axis


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
    obs_xaxis, obs_glo, obs_nh, obs_sh = generate_axes(obs, frac=args.frac)

    # model axis
    model = ori(args.inDir, args.beg_year, args.end_year, 'model')
    model_xaxis, model_glo, model_nh, model_sh = \
        generate_axes(model, frac=args.frac)

    # Configure the plot
    fig, (ax_glo, ax_nh, ax_sh) = \
        matplotlib.pyplot.subplots(3, 1,
                                   sharex=True,
                                   layout="constrained")
    fig.suptitle(f"{storm_type} Count By Latitude ({year_range})")

    # Global
    model_plt, = ax_glo.plot(model_xaxis, model_glo, label=args.expName)
    obs_plt, = ax_glo.plot(obs_xaxis, obs_glo, label="obs")
    ax_glo.set_title('Global', loc="left")
    ax_glo.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax_glo.set_ylim(0)
    ax_glo.axhline(y=0)
    ax_glo.grid(True, which="both")
    ax_glo.legend(handels=[model_plt, obs_plt])

    # Northern Hemisphere
    model_nh_plt, = ax_nh.plot(model_xaxis, model_nh, label=args.expName)
    obs_nh_plt, = ax_nh.plot(obs_xaxis, obs_nh, label="obs")
    ax_nh.set_title('Northern Hemisphere', loc="left")
    ax_nh.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax_nh.set_ylim(0)
    ax_nh.axhline(y=0)
    ax_nh.grid(True, which="both")
    ax_nh.set_ylabel(f"no. {storm_type} per year")

    # Southern Hemisphere
    model_sh_plt, = ax_sh.plot(model_xaxis, model_sh, label=args.expName)
    obs_sh_plt, = ax_sh.plot(obs_xaxis, obs_sh, label="obs")
    ax_sh.set_title('Southern Hemisphere', loc="left")
    ax_sh.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
    ax_sh.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(30))
    ax_sh.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(15))
    ax_sh.set_ylim(0)
    ax_sh.axhline(y=0)
    ax_sh.grid(True, which="both")
    ax_sh.set_ylabel(f"no. {storm_type} per year")
    ax_sh.set_xlabel(u"longitude (\N{DEGREE SIGN})")

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)

        plot_filename = "by_longitude.pdf"
        matplotlib.pyplot.savefig(plot_filename)
        shutil.copyfile(plot_filename,
                        os.path.join(args.outDir, plot_filename))
        print(f"Plot stored in '{os.path.join(args.outDir, plot_filename)}'")
