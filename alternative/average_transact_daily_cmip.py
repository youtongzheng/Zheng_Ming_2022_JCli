import numpy as np
import xarray as xr
import glob
import sys

# Parse inputs
if len(sys.argv) !=6:
    print('ERROR---Usage: python average_transact_daily_cmip.py ' +
        '<year> <hemisphere> <dir> <dir1>' )
    exit()
    
exp = sys.argv[1] 
year = sys.argv[2]
hemisphere = sys.argv[3] 
path = sys.argv[4]
path2 = sys.argv[5]

#
filenames =  glob.glob(path + hemisphere + "_*_daily_transect_" + exp +"_cmip.nc")
print(filenames)
f = xr.open_mfdataset(filenames,decode_cf = True)

filenames =  glob.glob(path2 + hemisphere + "_*_daily_transect_" + exp +".nc")
print(filenames)
ftmp = xr.open_mfdataset(filenames,decode_cf = True)

f = f.sel(time=slice("2006-06-12", "2014-12-31"))
ftmp = ftmp.sel(time=slice("2006-06-12", "2014-12-31"))

f['v10_ave'] = ftmp.v_ref.sel(distance = slice(-3.5, 3.5)).mean(dim = 'distance')
f['u10_ave'] = ftmp.u_ref.sel(distance = slice(-3.5, 3.5)).mean(dim = 'distance')

f = f.rename({'level': 'altitude'})

f = f.assign(temp=ftmp.temp)
f = f.assign(sphum=ftmp.sphum)

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

ffpos.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_onice.nc")
ffpos_weaku.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_onice_weaku.nc")
ffpos_strongu.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_onice_strongu.nc")

ffneg.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_office.nc")
ffneg_weaku.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_office_weaku.nc")
ffneg_strongu.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_office_strongu.nc")

ff.to_netcdf(path + hemisphere + '_' + year + "_average_transect_" + exp +"_cmip_all.nc")