import numpy as np
import xarray as xr
import glob
import sys
from num_utils import *

# Parse inputs
if len(sys.argv) !=5:
    print('ERROR---Usage: python Gen_datasets_distance_from_ice_edge_daily.py ' +
        '<year> <hemisphere> <dir>' )
    exit()
    
exp = sys.argv[1] 
year = sys.argv[2]
hemisphere = sys.argv[3] 
outpath = sys.argv[4]

if hemisphere == 'nh':
    lat0 = 70.
    lat1 = 85.
else:
    lat0 = -85.
    lat1 = -50. 

path = '/archive/Youtong.Zheng/awg/2021.02/' + exp + '/gfdl.ncrc3-intel21-prod-openmp/pp/atmos/ts/monthly/5yr/'

f = xr.open_dataset(path + 'atmos.' + year + '.qadt_dyn.nc')
f = f.isel(pfull = slice(17,33))

varnames3d = ['qadt_conv', 'qadt_vdif', 'qadt_lsform','qadt_super','qadt_eros', 'qadt_destr',
              'qldt_cond','qldt_liqadj','qldt_accr','qldt_auto','qldt_berg','qldt_conv','qldt_destr','qldt_dyn','qldt_eros',
              'qldt_evap','qldt_fill','qldt_freez','qldt_ls','qldt_rime','qldt_vdif',
              'tdt_lw', 'tdt_sw', 'qdt_conv', 'qdt_vdif',
              'cld_amt','ice_wat','liq_wat',
              'temp','sphum','omega', 'mc_full', 'z_full']

for i, var in enumerate(varnames3d):
    f0 = xr.open_dataset(path + 'atmos.' + year + '.' + var + '.nc')
    f0 = f0.isel(pfull = slice(17,33))
    f = f.assign({var:f0[var]})


varnames2d = ['low_cld_amt','lwdn_sfc_clr','lwdn_sfc', 'swdn_sfc_clr','swdn_sfc','z_pbl','evap']

for i, var in enumerate(varnames2d):
    f0 = xr.open_dataset(path + 'atmos.' + year + '.' + var + '.nc')
    f = f.assign({var:f0[var]})
    
path = "/archive/oar.gfdl.cmip6/CM4/warsaw_201803/CM4_amip/gfdl.ncrc4-intel16-prod-openmp/pp/atmos/ts/monthly/15yr/"
f0 = xr.open_mfdataset(glob.glob(path + '*ice_mask.nc'),
                      combine = 'by_coords',decode_cf = True)

#read ice
f0 = f0.sel(time=slice(year[0:4] + "-01-01", year[7:11] + "-12-31"))
f0 = f0.interp(lat=f.lat.values, lon=f.lon.values)
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
else:
    if hemisphere == 'sh_ross':
        f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-180, -130))
    else: 
        f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-50, 0))
        
# path = "/work/Youtong.Zheng/Daily/"
# tmp = xr.open_dataset(glob.glob(path + 'calipso/2D_phase/*20070101*.nc')[0])
# lat_plot = tmp.sel(latitude=slice(lat0, lat1)).latitude.values
# f0 = f0.interp(lat=lat_plot)

#extract values
lat = f0.lat.values
lon = f0.lon.values
time = f0.time.values 
lev = f0.pfull.values

nlat = len(lat)
nlon = len(lon)
ntime = len(time)

ice_achv = f0['ice_mask'].values
qadt_dyn_achv = f0['qadt_dyn'].values

print('start computing...')
qadt_dyn_out, ice_out = Binning_distance_from_ice_edge(ice_achv, qadt_dyn_achv, ntime, nlat, nlon, hemisphere, model = True)

ds_tmp = xr.DataArray(qadt_dyn_out, coords=[lev, np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level", "distance","time"]).rename('qadt_dyn')

ds_tmp1 = xr.DataArray(ice_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('ice_mask')

ds_out = xr.merge([ds_tmp, ds_tmp1])
    
for ivar, var in enumerate(varnames3d):
    var_achv = f0[var].values
    var_out, ice_out = Binning_distance_from_ice_edge(ice_achv, var_achv, ntime, nlat, nlon, hemisphere, model = True)
    
    ds_tmp = xr.DataArray(var_out, coords=[lev,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level","distance","time"]).rename(var)
    
    ds_out = xr.merge([ds_out, ds_tmp]) 

for ivar, var in enumerate(varnames2d):
    var_achv = f0[var].values
    var_out, ice_out = Binning_distance_from_ice_edge(ice_achv, var_achv, ntime, nlat, nlon, hemisphere, model = True)
    
    ds_tmp = xr.DataArray(var_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename(var)
    
    ds_out = xr.merge([ds_out, ds_tmp]) 
    
ds_out.attrs['title'] = 'year:' + year

# outpath = '/nbhome/Youtong.Zheng/'
ds_out.to_netcdf(outpath + hemisphere + '_' + year + "_monthly_transact_" + exp +".nc")
print('Job accomplished for year:' + year)
