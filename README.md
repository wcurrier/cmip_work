# cmip_work
Process and Analyze CMIP data for downscaling

Use cmip2icarForce.py to run subsetCMIP.py, which calls cmipFunctions.py

Example:

python -u /glade/u/home/currierw/cmip_work/cmip2icarForce.py -ins BCC -model bcc-csm1-1 -ensemble r1i1p1 -scen rcp45 > rcp45BCCout.txt
