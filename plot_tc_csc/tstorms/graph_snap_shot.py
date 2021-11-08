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

"""Tropical Storm Snapshot Plot Generator

This module generates the Tropical Storm Snapshot plot using tracking data
generated from a GCM.
"""

__all__ = [
    'generate_ori_data',
    'generate_plot_data',
    'write_plot_data',
]

import argparse
from .StormBox import read_storm_stats
from . import argparse as tsargparse
from .ori_stat import cat_ori_files
from .config import exeext, pkglibexecdir, gracebat
import os
import shutil
import subprocess
import jinja2
import tempfile


def generate_ori_data(type, template_env):
    ori_data = []
    with open('ori', 'r') as infile:
        for line in infile.readlines():
            ln = line.split()
            ori_data.append(f'{ln[0]} {ln[1]}')
    geog_dat = template_env.get_template('geog.dat')
    with open(f'ori_{type}.dat', 'w') as out:
        out.write(geog_dat.render(data=ori_data))


def generate_plot_data(statDir, oriDir, beg_year, end_year, type):
    # Namelist required for freq_ori command
    nml_input = b"""&input
do_40ns = .true.
do_map  = .false.
do_lon  = .true.
do_lat  = .false.
nexp    =  1
/
"""
    freq_ori_cmd = os.path.join(pkglibexecdir, 'freq_ori' + exeext)
    fstats = read_storm_stats(os.path.join(statDir, f'stats_{type}_{beg_year}-{end_year}'))

    xts = {}
    for b in ['NH', 'SH']:
        xts[b] = []
        for y in fstats[b].years:
            xts[b].append(f"{y} {fstats[b].get_year_total(y)}")

    # mean data
    xscyc = {}
    # Northern Hemisphere
    xscyc['NH'] = []
    for i, v in enumerate(fstats['NH'].mean[:12], start=1):
        xscyc['NH'].append(f"{i} {v}")
    # Souther Hemisphere
    xscyc['SH'] = []
    for i, v in enumerate(fstats['SH'].mean[6:12] + fstats['SH'].mean[:6], start=1):
        xscyc['SH'].append(f"{i} {v}")

    cat_ori_files(oriDir, beg_year, end_year)
    write_plot_data("xts_nh.dat", xts['NH'])
    write_plot_data("xts_sh.dat", xts['SH'])
    write_plot_data("xscyc_nh.dat", xscyc['NH'])
    write_plot_data("xscyc_sh.dat", xscyc['SH'])

    subprocess.run([freq_ori_cmd], input=nml_input)
    for region in ['gl', 'nh', 'sh']:
        with open(f"xlon_{region}.dat", 'a') as outfile:
            with open(f"flon_{region}") as infile:
                outfile.write(infile.read())

    return fstats['NH'].mean[12], fstats['SH'].mean[12]


def write_plot_data(file, array):
    with open(file, "a") as f:
        f.write("\n".join(array + ['&\n']))


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-o",
                           help="Directory where plots will be stored",
                           metavar="outDir",
                           dest="outDir",
                           default=os.getcwd(),
                           type=tsargparse.absPath,
                           action=tsargparse.createDir)
    argparser.add_argument("inDir",
                           help="Directory where tropical storm data are available",
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

    grace_template_dir = os.path.join(os.path.dirname(__file__), 'templates')

    template_env = jinja2.Environment(loader=jinja2.FileSystemLoader(grace_template_dir))
    template_env.keep_trailing_newline = True
    template_env.trim_blocks = True
    template_env.lstrip_blocks = True
    template_env.rstrip_blocks = True

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        nh_obs_mean, sh_obs_mean = generate_plot_data(args.inDir, args.obsDir, args.beg_year, args.end_year, 'obs')
        generate_ori_data('obs', template_env)
        nh_model_mean, sh_model_mean = generate_plot_data(args.inDir, args.inDir, args.beg_year, args.end_year, 'model')
        generate_ori_data('model', template_env)

        snap_shot_par = template_env.get_template('snap_shot.par')
        snap_shot_data = {
            "BEG_YEAR": args.beg_year,
            "END_YEAR": args.end_year,
            "NH_OBS_MEAN": nh_obs_mean,
            "NH_MODEL_MEAN": nh_model_mean,
            "SH_OBS_MEAN": sh_obs_mean,
            "SH_MODEL_MEAN": sh_model_mean,
            "PLOT_TITLE": args.expName,
        }
        with open('snap_shot.par', 'w') as out:
            out.write(snap_shot_par.render(snap_shot_data))

        grace_cmd = [
            gracebat,
            "-printfile", "snapshot.ps",
            "-param", "snap_shot.par",
            "-hardcopy",
            "-graph", "3", "xts_nh.dat", "-graph", "7", "xts_sh.dat",
            "-graph", "2", "xscyc_nh.dat", "-graph", "6", "xscyc_sh.dat",
            "-graph", "1", "xlon_nh.dat", "-graph", "5", "xlon_sh.dat",
            "-graph", "0", "ori_obs.dat", "-graph", "4", "ori_model.dat",
        ]
        subprocess.run(grace_cmd)

        shutil.copyfile('snapshot.ps',
                        os.path.join(args.outDir,
                                     'snap_shot_{0:04d}-{1:04d}.ps'.format(args.beg_year,
                                                                           args.end_year)))