#!/usr/bin/env python
# coding: utf-8

"""
SYNOPSIS

    cmip2icarForce.py [-h] [-ins [Institute]] [-model [MODEL]] [-ensemble [ensemble]]
                      [-scen [SCEN]] [-v] [--verbose]

DESCRIPTION

    Script to download data CMIP5 forcing data for ICAR by selecting model and time period.
    Specify GCM, time period, ensemble memeber, latitutde and longitude, working directory
    Sea Surface Temperature are regridded from rotated pole if necessary - fills from daily to 6 hourly
    U,V data are regridded to to a non-staggered grid (ACCESS and HadGEM)
    Converts longitude from 0-360 to -180-180
    Calculate model level heights or pressure levels depending on staggered
    Computes water vapor micing ratio from specific humidity

    TODO: Add convective precipitation

EXAMPLES

    cmip2icarForce.py ins BCC -model bcc-csm-1 -ensemble r1i1p1 -scen historical

EXIT STATUS

    TODO: List exit codes

AUTHOR

    William Currier - currierw@ucar.edu - wrote progrm to process CMIP data.
    Ethan Gutmann   - gutmann@ucar.edu  - wrote program to download CMIP data
                                          from CEDA

LICENSE

    This script is in the public domain.

VERSION
    1.0

"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse
import pandas as pd

global verbose
verbose=False

sys.path.append('/glade/u/home/currierw/cmip_work')
import subsetCmip

##########################################################################################################################
############################## CHANGE THE CODE BELOW IF CHANGING THINGS ##################################################
##########################################################################################################################



##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
############################## CHANGE THE CODE ABOVE IF CHANGING THINGS ##################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
####################################### ADDITIONAL INFORMATION ###########################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#
#                   [{"var_name":"va", "domain":"atmos", "interval":"6hr"},
#                   {"var_name":"ua", "domain":"atmos", "interval":"6hr"},
#                   {"var_name":"ta", "domain":"atmos", "interval":"6hr"},
#                   {"var_name":"ps", "domain":"atmos", "interval":"6hr"},
#                   {"var_name":"hus", "domain":"atmos", "interval":"6hr"},
#                   {"var_name":"prc", "domain":"atmos", "interval":"3hr"},
#                   {"var_name":"tos", "domain":"ocean", "interval":"day"},
#                   {"var_name":"orog", "domain":"atmos", "interval":"fx"},
#                   {"var_name":"sftlf", "domain":"atmos", "interval":"fx"},
#
# DESIRED MODEL LIST
# GCMList=['ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','BNU-ESM','CanESM2','CCSM4','CESM1-BGC', \
#     'CESM1-CAM5','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-g2','FIO-ESM', \
#     'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H-cc','GISS-E2-H','GISS-E2-R-cc','GISS-E2-R', \
#     'GISS-E2-R','HadGEM2-AO','HadGEM2-CC','HadGEM2-ES','HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-LR', \
#     'IPSL-CM5A-LR','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM','MIROC-ESM-CHEM','MIROC5', \
#     'MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']
#
# MODEL LIST ON CHEYENNE
# chyModelList=['BCC/bcc-csm1-1','BCC/bcc-csm1-1-m','BNU/BNU-ESM','CCCma/CanESM2','CNRM-CERFACS/CNRM-CM5','CSIRO-BOM/ACCESS1-0','CSIRO-BOM/ACCESS1-3', \
#               'CSIRO-QCCCE/CSIRO-Mk3-6-0','IPSL/IPSL-CM5A-LR','IPSL/IPSL-CM5A-MR','IPSL/IPSL-CM5B-LR', \
#               'LASG-CESS/FGOALS-g2','MIROC/MIROC5','MIROC/MIROC-ESM-CHEM','MIROC/MIROC-ESM','MOHC/HadGEM2-ES','MRI/MRI-CGCM3',
##########################################################################################################################
################################################## OTHER NOTES ###########################################################
##########################################################################################################################
#
#               'NCC/NorESM1-M','NOAA-GFDL/GFDL-CM3','NOAA-GFDL/GFDL-ESM2G']
# If it is possible to get all of the MPI data (with divergence / vorticity winds) these can be converted
#  to u and v winds using something like "cdo -s -r -f nc setreftime,1860-01-01,00:00 -settunits,day -sp2gp -dv2uv input.grb output.nc"
#
#  to remap a rotated pole (e.g. ocean temperatures) first create a cdo grid description:
#   cdo griddes sample_file.nc  > grid.txt
#  then remap to it
#   cdo -f nc -s remapcon,grid.txt rotated_pole_data.nc output_data.nc
#
# ftp://ftp.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/6hrLev/ps
#
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

# def main(ins=ins,model=model,ensemble=ensemble,scen=scen):
#         try:
#             subsetCmip.process(model=model,ins=ins,ensemble=ensemble,scen=scen)
#         except Exception as e:
#             print("Error: ")
#             print(e)

# if __name__ == '__main__':
    # main()
try:
    # initiate the parser
    parser= argparse.ArgumentParser(description='Download ICAR data from CMIP5 FTP archive. ',
                                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-mip',   dest="mip",       nargs="?", action='store', default="cmip5",       help="MIP to look for data from (cmip5,cmip6)")
    parser.add_argument('-ins',   dest="ins",       nargs="?", action='store', default="NCAR",     help="CMIP5 Model Institution (e.g. BCC, NCAR, ...)")
    parser.add_argument('-model', dest="model",     nargs="?", action='store', default="CCSM4",     help="CMIP5 Model name (e.g. bcc-csm1-1, CCSM4, ...)")
    parser.add_argument('-ensemble',   dest="ensemble",       nargs="?", action='store', default="r6i1p1",      help="run, initialization, physics member (e.g. r1i1p1, r6i1p1)")
    parser.add_argument('-scen',  dest="scen",      nargs="?", action='store', default="historical",  help="scenario (e.g. historical, rcp45, rcp85)")
    parser.add_argument('-v', '--version',action='version',
                            version='CMIP5 ICAR ingest 1.0')
    parser.add_argument ('--verbose', action='store_true',
                            default=False, help='verbose output', dest='verbose')
    args = parser.parse_args()

    ins=args.ins
    print(ins)
    model=args.model
    ensemble=args.ensemble
    scen=args.scen

    # exit_code = main(ins=ins,model=model,ensemble=ensemble,scen=scen)
    subsetCmip.process(model=model,ins=ins,ensemble=ensemble,scen=scen)

    if exit_code is None:
        exit_code = 0
        sys.exit(exit_code)
except KeyboardInterrupt as e: # Ctrl-C
    raise e
except SystemExit as e: # sys.exit()
    raise e
except Exception as e:
    print('ERROR, UNEXPECTED EXCEPTION')
    print(str(e))
    traceback.print_exc()
    os._exit(1)
