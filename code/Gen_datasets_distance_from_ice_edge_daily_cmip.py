import numpy as np
import xarray as xr
import glob
import sys
from num_utils import *

# Parse inputs
if len(sys.argv) !=5:
    print('ERROR---Usage: python Gen_datasets_distance_from_ice_edge_daily_cmip.py ' +
        '<year> <hemisphere> <dir>' )
    exit()
    
exp = sys.argv[1] 
year = sys.argv[2]
hemisphere = sys.argv[3] 
outpath = sys.argv[4]

if hemisphere[0:2] == 'nh':
    lat0 = 70.
    lat1 = 85.
else:
    lat0 = -85.
    lat1 = -50. 

#COSP data
path = '/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos_cmip_2deg_daily_COSP/ts/daily/6yr/'

f = xr.open_dataset(path + 'atmos_cmip_2deg_daily_COSP.' + year + '.clcalipso.nc')

#CMIP 2D data
path = '/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos_cmip_2deg_daily_2D/ts/daily/6yr/'

varnames2d = ['cllcalipso','vas', 'uas', 'ts']
for i, var in enumerate(varnames2d):
    f0 = xr.open_dataset(path + 'atmos_cmip_2deg_daily_2D.' + year + '.' + var + '.nc')
    f = f.assign({var:f0[var]})
    
#CMIP 3D data
path = '/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos_cmip_2deg_daily_3D/ts/daily/6yr/'

varnames3d = ['ta','hus','hur','wap']
for i, var in enumerate(varnames3d):
    f0 = xr.open_dataset(path + 'atmos_cmip_2deg_daily_3D.' + year + '.' + var + '.nc')
    f = f.assign({var:f0[var]})

#ice data    
path = "/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos/ts/monthly/15yr/"
f0 = xr.open_mfdataset(glob.glob(path + '*ice_mask.nc'),
                      combine = 'by_coords',decode_cf = True)

#read ice
f0 = f0.sel(time=slice(year[0:4] + "-01-01", year[9:13] + "-12-31"))
f0 = f0.interp(lat=f.lat.values, lon=f.lon.values)

f0 = f0.reindex_like(f.time.sel(time=f.time.dt.day == 1), method = 'bfill')
f0 = f0.reindex_like(f.time, method = 'ffill')
f0['time'] = f.time
f = f.assign(ice_mask=f0.ice_mask)

#read land mask
path = "/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos/av/annual_15yr/"
f0 = xr.open_mfdataset(glob.glob(path + '*.nc'),
                      combine = 'by_coords',decode_cf = True)
f0 = f0.interp(lat=f.lat.values, lon=f.lon.values)
f = f.where(f0.land_mask < 0.5)

#pre-process
f.coords['lon'] = (f.coords['lon'] + 180) % 360 - 180
f = f.sortby(f.lon)

if hemisphere == 'nh':
    f0 = f.sel(lat=slice(lat0, lat1)).sel(lon=slice(-15, 55))
elif hemisphere == 'nh-chukchi':
    f0 = xr.concat([f.sel(lat=slice(lat0, lat1)).sel(lon=slice(-180, -165)),
                    f.sel(lat=slice(lat0, lat1)).sel(lon=slice(170, 180))], dim="lon")
elif hemisphere == 'sh':
    f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-50, 0))
else:
    f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-180, -130))

#interpolate to the calipso grids to make sure that the influence of different horizontal grids is minimized       
path = "/work/Youtong.Zheng/Data/satellite/"
tmp = xr.open_dataset(glob.glob(path + 'calipso/2D_phase/*20070101*.nc')[0])
lat_plot = tmp.sel(latitude=slice(lat0, lat1)).latitude.values
f0 = f0.interp(lat=lat_plot, lon=f.lon.values)

#extract values
lat = f0.lat.values
lon = f0.lon.values
time = f0.time.values 
altitude = f0.alt40.values/1000.
plev19 = f0.plev19.values/100.

nlat = len(lat)
nlon = len(lon)
ntime = len(time)

ice_achv = f0['ice_mask'].values
clcalipso_achv = f0['clcalipso'].values

print('start computing...')
clcalipso_out, ice_out = Binning_distance_from_ice_edge(ice_achv, clcalipso_achv, ntime, nlat, nlon, hemisphere[0:2], model = True)

ds_tmp = xr.DataArray(clcalipso_out, coords=[altitude, np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["altitude", "distance","time"]).rename('clcalipso')

ds_tmp1 = xr.DataArray(ice_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('ice_mask')

ds_out = xr.merge([ds_tmp, ds_tmp1])
    
for ivar, var in enumerate(varnames3d):
    var_achv = f0[var].values
    var_out, ice_out = Binning_distance_from_ice_edge(ice_achv, var_achv, ntime, nlat, nlon, hemisphere, model = True)
    
    ds_tmp = xr.DataArray(var_out, coords=[plev19,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["plev19","distance","time"]).rename(var)
    
    ds_out = xr.merge([ds_out, ds_tmp]) 

for ivar, var in enumerate(varnames2d):
    var_achv = f0[var].values
    var_out, ice_out = Binning_distance_from_ice_edge(ice_achv, var_achv, ntime, nlat, nlon, hemisphere[0:2], model = True)
    
    ds_tmp = xr.DataArray(var_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename(var)
    
    ds_out = xr.merge([ds_out, ds_tmp]) 
    
ds_out.attrs['title'] = 'year:' + year

ds_out.to_netcdf(outpath + hemisphere + '_' + year + "_daily_transect_" + exp +"_cmip.nc")
print('Job accomplished for year:' + year)
