# cmip_work

## Use: 

   Script to download and subset CMIP5 forcing data for ICAR
   Forcing data include 6 hour - 3-Dimensional data:
   * U,V wind speeds
   * Temperature
   * Specific Humidity
   * Surface Pressure (2D)
   * Sea Surface Temperature (2D - Rotated Pole - daily -> 6 hour -> interpolated)
   * Orography (2D - fixed in time)
   Select model and GCM experiment (historical/rcp45/rcp85).
   Specify institute, GCM, experiment, ensemble memeber. 
 
 ## Processing:
 
   * Sea Surface Temperature are regridded from rotated pole if necessary - fills from daily to 6 hourly - uses CDO commands. See `transformTOS.ipynb` for example.
   * Sea Surface Temperature are interpolated longitudnally to fill misising values over land.
   * U,V data are regridded to to a non-staggered grid (ACCESS and HadGEM modesl)
   * Converts longitude from 0-360 to -180-180
   * Calculate model level heights or pressure levels depending on GCM. ICAR v2.0 can read either and compute the rest
   * Computes water vapor mixing ratio from specific humidity
 
 ## Notes:
 
   * Currently set up for personal storage on NCAR's Cheyenne computer.
   * Writes out intermediate, spatially subsetted, 10 year files - stitches them back together
   * Latitude and longitude are hadcoded at the top of subsetCMIP.py - currently set up for western US.
   * Will process entire historic, rcp45, and rcp85 periods depending on specification - hardcoded for now to make inputs smaller
   * Need `get_dataset` function from https://github.com/gutmann/cmip_ingest/tree/master/scripts/download 
   * Accesses CEDA website: ftp://ftp.ceda.ac.uk/badc/ need to register an acccount with them. Good idea to create a .netrc file with username and password.
   * Using NCAR's Cheyenne computer? See what Cheyenne has for 6 hour model level data  with `exploreCMIPonCheyenne.ipynb`. Note that only historical 6 hour model level data exists on Cheyenne. 
    
## TODO: 
    
   Add convective precipitation
   Add smarter way to destagger ACCESS and HADGEM - currently interpolates grid using xr.interp

## Example:

Use cmip2icarForce.py to run subsetCMIP.py, which calls cmipFunctions.py

    python -u /glade/u/home/currierw/cmip_work/cmip2icarForce.py -ins BCC -model bcc-csm1-1 -ensemble r1i1p1 -scen rcp45 > rcp45BCCout.txt
