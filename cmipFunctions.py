#!/usr/bin/env python
# coding: utf-8
# In[ ]:

# cmipFunctions.py

# Functions needed to convert cmip5 data to subsetted and merged files to run ICAR

# convertTOS
#     Convert Sea Surface Temperature files (TOS) from the rotated pole grid to
#     regular grid that all other variables are on using cdo commmands:
#     https://code.mpimet.mpg.de/projects/cdo/embedded/index.html

# checkFiles
#     Check that files exist for variables and that we have enough files
#     to be able to put together a dataset to run ICAR over 150 year period

# loadNonStaggeredVars
#     Loads in, subsets spatailly and temporally, writes out subsetted files
#     Only files that are not should be used. Staggered varaibles should be further
#     and interpolated in x and y.


# Created by William Ryan Currier on 2020-01-27.
import numpy as np
import pandas as pd
import xarray as xr
import glob
import sys
from cdo import *
import os
cdo = Cdo()


def convertTos(outDir,inDir):
    # Check to see if files exist, if they don't let's regrid them using cdo commands
    regrdFiles=glob.glob(outDir+'tos*regrd.nc')
    if len(regrdFiles) == 0 :
        # Get the grid dimensions, lat, and lon of the grid we want to transform to (tmp.nc)
        cdo.griddes("-f "+outDir+'tmp.nc'+" >"+outDir+"grid.txt")

        fNames=sorted(glob.glob(inDir+'tos*.nc')) # get all base sea surface temperature variables
        for i , f in enumerate(fNames):
            outputBase=os.path.splitext(fNames[i].split('/')[-1])[0] # removes path and file extension of each file
            cdo.remapbil("grid.txt ", input=f, output=outDir+outputBase+"_regrd.nc", options="-f nc")
            print("Created new TOS file "+outDir+outputBase+"_regrd.nc")

def checkFiles(inDir,variable,time_period,chyColl):

    # Cheyenne Collection has each variable stored in it's own directory
    if chyColl:
        inDir=inDir+variable
        dirVar=inDir+'/'+variable+'*'
    else:
        dirVar=inDir+variable+'*'

    # Get time period difference from filename
    if variable == 'tos':
        file1 = sorted(glob.glob(dirVar+'*regrd.nc'))
    #         fileParts=file1[0].split("_")     # underscore before dates
    #         yearStrs=fileParts[-2].split("-") # break apart years from dates, which is the last element in file string
    #         # Convert to date string
    #         yearStart1=pd.to_datetime(yearStrs[0],format='%Y%m%d')
    #         yearStart2=pd.to_datetime(yearStrs[1],format='%Y%m%d')
    else:
        file1 = sorted(glob.glob(dirVar+'*.nc'))
    #         fileParts=file1[0].split("_")     # underscore before dates
    #         yearStrs=fileParts[-1].split("-") # break apart years from dates, which is the last element in file string
    #         yearStr1=yearStrs[0]              # first year
    #         yearStrs2=yearStrs[1].split(".")  # remove .nc extension
    #         yearStr2=yearStrs2[0]             # second year - remove .nc extension
    #         # Convert to date string
    #         yearStart1=pd.to_datetime(yearStr1,format='%Y%m%d%H')
    #         yearStart2=pd.to_datetime(yearStr2,format='%Y%m%d%H')

    #     fTPeriod = (yearStart2-yearStart1)

    #     oneYearDelt = pd.to_datetime('2006010400',format='%Y%m%d%H')-pd.to_datetime('2005010100',format='%Y%m%d%H') # make longer than 365, leap year

    # Good for variables Hus, va, ua, ta
    if variable == 'hus' or variable == 'va' or variable == 'ua' or variable == 'ta':
        if time_period == 'historical':
            if len(file1)<56 and chyColl:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)
        if time_period == 'rcp45' and chyColl:
            if len(file1)<94:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)
        if time_period == 'rcp85' and chyColl:
            if len(file1)<94:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)
    if variable == 'ps' or variable == 'tos':
        if time_period == 'historical' and chyColl:
            if len(file1)<2:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)
        if time_period == 'rcp45' and chyColl:
            if len(file1)<2:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)
        if time_period == 'rcp85' and chyColl:
            if len(file1)<2:
                sys.exit("Likely need more data files for "+variable)
            else:
                print("All files available for "+variable)

def loadNonStaggeredVarsResample(inDir,outDir,variable,Model,time_period,ensemble,lat_bnds,lon_bnds,startDates,endDates,chyColl,rsTime,rsMet):
    # 1. Load in rectilinear data
    # 2. Subset the data spatially
    # 3. Subset by time
    # 4. Write out hourly subsetted file in 10 year periods

    # inDir       = directory where datafiles are located - string
    # outDir      = where to write the datafiles
    # variable    = variable we want to subset - string
    # Model       = What GCM? - string
    # time_period = historical,rcp45,rcp85 - string
    # ensemble    = r1p1i1 - string
    # lat_bnds    = latitutde boundaries [south,north] -tuple
    # lon_bnds    = longitude boundaries [west,east] - tuple
    # startDates  = list of start times to slice the data into - datetime
    # endDates    = list of end dates to slice the data into - datetime
    # chyColl     = cheyenne collection [True/False] - do dat already exist there?
    # rsTime      = time you want to resample to - '6hr'
    # rsMet       = resampling method - e.g. 'sum' - need to add more functionality if needed


    if chyColl:
        inDirUp=inDir+variable+'/'
    else:
        inDirUp=inDir

    # Load in rectilinear data
    ds = xr.open_mfdataset(inDirUp+variable+'*.nc',combine='by_coords',engine='netcdf4',parallel=True) # load data
    # Subset the data spatially
    ds = ds.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
    # Resample
    if rsMet=='sum' and rsTime == '6hr':
        ds  = ds['prc'].resample(time='6H').sum()
    # Subset by time and write out hourly file in 10 year periods
    for startDateLs, endDateLs in zip(startDates, endDates):
        ds_sub_Time  = ds.sel(time=slice(startDateLs.strftime('%Y-%m-%dT%H:%M:%S'),endDateLs.strftime('%Y-%m-%dT%H:%M:%S')))  # Time slice
        ds_sub_Time.to_netcdf(outDir+variable+'_6hrLev_' + Model +'_'+ time_period +'_'+ ensemble +'_' \
                                + startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d') \
                                +'_'+ ensemble + '_subset.nc')
        print("Writing out "+outDir+variable+'_6hr_' + Model +'_'+startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d')+'_'+ ensemble+ '_subset.nc')
    del ds
    del ds_sub_Time

def loadNonStaggeredVars(inDir,outDir,variable,Model,time_period,ensemble,lat_bnds,lon_bnds,startDates,endDates,chyColl):
    # 1. Load in rectilinear data
    # 2. Subset the data spatially
    # 3. Subset by time
    # 4. Write out hourly subsetted file in 10 year periods

    # inDir       = directory where datafiles are located - string
    # outDir      = where to write the datafiles
    # variable    = variable we want to subset - string
    # Model       = What GCM? - string
    # time_period = historical,rcp45,rcp85 - string
    # ensemble    = r1p1i1 - string
    # lat_bnds    = latitutde boundaries [south,north] -tuple
    # lon_bnds    = longitude boundaries [west,east] - tuple
    # startDates  = list of start times to slice the data into - datetime
    # endDates    = list of end dates to slice the data into - datetime
    # chyColl     = cheyenne collection [True/False] - do dat already exist there?
    if chyColl:
        inDirUp=inDir+variable+'/'
    else:
        inDirUp=inDir

    # Load in rectilinear data
    ds = xr.open_mfdataset(inDirUp+variable+'*.nc',combine='by_coords',engine='netcdf4',parallel=True) # load data
    # Subset the data spatially
    ds = ds.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
    # Subset by time and write out hourly file in 10 year periods
    for startDateLs, endDateLs in zip(startDates, endDates):
        ds_sub_Time  = ds.sel(time=slice(startDateLs.strftime('%Y-%m-%dT%H:%M:%S'),endDateLs.strftime('%Y-%m-%dT%H:%M:%S')))  # Time slice
        ds_sub_Time.to_netcdf(outDir+variable+'_6hrLev_' + Model +'_'+ time_period +'_'+ ensemble +'_' \
                                + startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d') \
                                +'_'+ ensemble+ '_subset.nc')
        print("Writing out "+outDir+variable+'_6hrLev_' + Model +'_'+startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d')+'_'+ ensemble+ '_subset.nc')
    del ds
    del ds_sub_Time


def loadStaggeredVars(inDir,outDir,variable,Model,time_period,ensemble,lat_bnds,lon_bnds,startDates,endDates,chyColl):
    # 1. Load in rectilinear - staggered data
    # 2. Alter the latitude and longitude boundaries
    # 3. Subset the data spatially
    # 4. Interpolate the data to the ta grid
    # 5. Subset by time
    # 6. Write out hourly subsetted file in 10 year periods

    # inDir       = directory where datafiles are located - string
    # outDir      = where to write the datafiles
    # variable    = variable we want to subset - string
    # Model       = What GCM? - string
    # time_period = historical,rcp45,rcp85 - string
    # ensemble    = r1p1i1 - string
    # lat_bnds    = latitutde boundaries [south,north] -tuple - same as for non-staggered data
    # lon_bnds    = longitude boundaries [west,east] - tuple - same as for non-staggered data
    # startDates  = list of start times to slice the data into - datetime
    # endDates    = list of end dates to slice the data into - datetime
    # chyColl     = cheyenne collection [True/False] - do dat already exist there?
    if chyColl:
        inDirUp=inDir+variable+'/'
    else:
        inDirUp=inDir

    # Load in rectilinear - possibly staggered data (u,v winds)
    ds = xr.open_mfdataset(inDirUp+variable+'*.nc',combine='by_coords',engine='netcdf4',parallel=True) # load data

    # Load Temperature File - subsetted
    if chyColl:
        tFiles=sorted(glob.glob(inDir+'ta/ta*.nc'))
    else:
        tFiles=sorted(glob.glob(inDir+'ta*.nc'))
        tFilesSub=sorted(glob.glob(inDir+'ta*subset.nc'))
        if len(tFilesSub)>0:
            for f in tFilesSub:
                tFiles.remove(f)
    dsT=xr.open_dataset(tFiles[0])

    # Check to see if variables are staggered
    if variable == 'ua':
        if np.abs(ds['lon'][0]-dsT['lon'][0])>0:
            staggered=True
        else:
            staggered=False
    elif variable == 'va':
        if np.abs(ds['lat'][0]-dsT['lat'][0])>0:
            staggered=True
        else:
            staggered=False
    if staggered:
        print("Data are staggered")
        # Alter the latitude and longitude boundaries
        dx=(ds.lon[1]-ds.lon[0])/2
        dy=(ds.lat[1]-ds.lat[0])/2
        # Slciing for V winds (north-south) add dy degree buffer to south and north
        lat_bndsV,lon_bndsV = [lat_bnds[0]-dy,lat_bnds[1]+dy], [lon_bnds[0], lon_bnds[1]]
        # Slicing for U winds (east-west) add dx degree buffer to west and east
        lat_bndsU,lon_bndsU = [lat_bnds[0],lat_bnds[1]],     [lon_bnds[0]-dx, lon_bnds[1]+dx]
        # Load in Temperature Data to interpolate - subsetted temperature
        tFilesSub=sorted(glob.glob(outDir+'ta*.nc'))
        dsTSub=xr.open_dataset(tFilesSub[0])
        # Subset and Interpolate the file
        if variable == 'ua':
            ds = ds.sel(lat=slice(*lat_bndsU), lon=slice(*lon_bndsU))
            if len(ds['lon'])-len(dsTSub['lon'])==1:
                ds = ds.interp(lat=dsTSub['lat'], lon=dsTSub['lon'])
            else:
                sys.exit('U interpolation needs inspection')
        elif variable == 'va':
            ds = ds.sel(lat=slice(*lat_bndsV), lon=slice(*lon_bndsV))
            if len(ds['lat'])-len(dsTSub['lat'])==1:
                ds = ds.interp(lat=dsTSub['lat'], lon=dsTSub['lon'])
            else:
                sys.exit('V interpolation needs inspection')
    else:
        print("Data are not staggered")
        # Subset Spaitally
        ds = ds.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))

    # Subset by time and write out hourly file in 10 year periods
    for startDateLs, endDateLs in zip(startDates, endDates):
        ds_sub_Time  = ds.sel(time=slice(startDateLs.strftime('%Y-%m-%dT%H:%M:%S'),endDateLs.strftime('%Y-%m-%dT%H:%M:%S')))  # Time slice
        ds_sub_Time.to_netcdf(outDir+variable+'_6hrLev_' + Model +'_'+ time_period +'_'+ ensemble +'_' \
                                + startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d') \
                                +'_'+ ensemble+ '_subset.nc')
        print("Writing out "+outDir+variable+'_6hrLev_' + Model +'_'+startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d')+'_'+ ensemble+ '_subset.nc')
    del ds
    del ds_sub_Time


def processTOS(inDir,outDir,variable,Model,time_period,ensemble,lat_bnds,lon_bnds,startDates,endDates):
    # 1. Load in regrided rectilinear TOS data
    # 2. Linearly interpolate ocean temperature - longitundally
    # 3. Subset spatially
    # 4. Resample to 6 hourly data - ffill
    # 5. Write out hourly subsetted file in 10 year periods

    # Load in data
    ds = xr.open_mfdataset(inDir+variable+'*_regrd.nc',combine='by_coords',engine='netcdf4',parallel=True) # load data
    # Interpolate to fill in NaN values - do this longitudnally - dont use subset as it wont fill in land not between oceans/seas
    ds['tos'] = ds['tos'].interpolate_na('lon',method='linear')
    # Subset Spaitally
    ds_sub = ds.sel(lat=slice(*lat_bnds), lon=slice(*lon_bnds))
    del ds
    # Resample Sea Surface Temp Data to 6 Hour data from Daily Mean - does forward fill
    ds6Hour_sub  = ds_sub.resample(time='6H').ffill()
    del ds_sub
    # Subset by time and write out hourly file in 10 year periods
    for startDateLs, endDateLs in zip(startDates, endDates):
        ds_sub_Time  = ds6Hour_sub.sel(time=slice(startDateLs.strftime('%Y-%m-%dT%H:%M:%S'),endDateLs.strftime('%Y-%m-%dT%H:%M:%S')))  # Time slice
        ds_sub_Time.to_netcdf(outDir+variable+'_6hr_' + Model +'_'+ time_period +'_'+ ensemble +'_'\
                                + startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d') \
                                +'_'+ ensemble+ '_regrd_subset.nc')
        print("Writing out "+outDir+variable+'_6hr_' + Model +'_'+ time_period +'_'+ ensemble +'_'\
                                + startDateLs.strftime('%Y%m%d')+'-'+endDateLs.strftime('%Y%m%d') \
                                +'_'+ ensemble+ '_regrd_subset.nc')
        del ds_sub_Time
    del ds6Hour_sub


def searchChyColl(fName,Model,variable):
    # find /glade/collections/cmip/cmip5/output1/*/*/*/6hr/atmos/6hrLev/r1i1p1/*/* -type d > fName
    # fName     = file with list of directories of 6 hour atmospheric model level data
    # Model     = string of model name e.g 'bcc-csm-1'
    # variable  = string with model variable name, e.g. 'hus'

    fList=[]                     # List of all files
    search = open(fName)
    for line in search:
        if Model+'/' in line:
            line = line.rstrip() # remove trailing \n
            fList.append(line)

    varDir=[]
    if len(fList)>0:
        for line2 in fList:
            if variable in line2:
                varDir.append(line2)

    if len(varDir)>0:
        files=sorted(glob.glob(str(varDir[0]+'/*.nc')))
    else:
        files=[]

    return fList, varDir, files

