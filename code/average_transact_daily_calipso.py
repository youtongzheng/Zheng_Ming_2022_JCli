import numpy as np
import xarray as xr
import glob
import sys

# Parse inputs
if len(sys.argv) !=5:
    print('ERROR---Usage: python average_transact_daily_calipso.py ' +
        '<year> <hemisphere> <dir>' )
    exit()
    
exp = sys.argv[1] 
year = sys.argv[2]
hemisphere = sys.argv[3] 
path = sys.argv[4]

f = xr.open_mfdataset(glob.glob(path + hemisphere + '_*_daily_transect_' + exp + '.nc'),
                                decode_cf = True)

f['v10_ave'] = f.v10.sel(distance = slice(-3.5, 3.5)).mean(dim = 'distance')
f['u10_ave'] = f.u10.sel(distance = slice(-3.5, 3.5)).mean(dim = 'distance')

if hemisphere[0:2] == 'nh':
    fpos = f.where(f.v10_ave >= 3, drop = True)
    fpos_weaku = f.where((f.v10_ave >= 3) & (f.u10_ave >= -2) & (f.u10_ave <= 2), drop = True)
    fpos_strongu = f.where((f.v10_ave >= 3) & ((f.u10_ave < -2) | (f.u10_ave > 2)), drop = True)

    fneg = f.where(f.v10_ave < -3, drop = True)
    fneg_weaku = f.where((f.v10_ave < -3) & (f.u10_ave >= -2) & (f.u10_ave <= 2), drop = True)
    fneg_strongu = f.where((f.v10_ave < -3) & ((f.u10_ave < -2) | (f.u10_ave > 2)), drop = True)
else:
    fpos = f.where(f.v10_ave <= -3, drop = True)
    fpos_weaku = f.where((f.v10_ave <= -3) & (f.u10_ave >= -2) & (f.u10_ave <= 2), drop = True)
    fpos_strongu = f.where((f.v10_ave <= -3) & ((f.u10_ave < -2) | (f.u10_ave > 2)), drop = True)

    fneg = f.where(f.v10_ave > 3, drop = True)
    fneg_weaku = f.where((f.v10_ave > 3) & (f.u10_ave >= -2) & (f.u10_ave <= 2), drop = True)
    fneg_strongu = f.where((f.v10_ave > 3) & ((f.u10_ave < -2) | (f.u10_ave > 2)), drop = True) 
    
ffpos = fpos.groupby("time.season").mean()
ffpos_weaku = fpos_weaku.groupby("time.season").mean()
ffpos_strongu = fpos_strongu.groupby("time.season").mean()

ffneg = fneg.groupby("time.season").mean()
ffneg_weaku = fneg_weaku.groupby("time.season").mean()
ffneg_strongu = fneg_strongu.groupby("time.season").mean()

ff = f.groupby("time.season").mean()

ffpos.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_onice.nc")
ffpos_weaku.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_onice_weaku.nc")
ffpos_strongu.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_onice_strongu.nc")

ffneg.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_office.nc")
ffneg_weaku.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_office_weaku.nc")
ffneg_strongu.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_office_strongu.nc")

ff.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_all.nc")