import numpy as np
import xarray as xr
import glob
import sys
from num_utils import *

# Parse inputs
if len(sys.argv) !=4:
    print('ERROR---Usage: python Gen_datasets_distance_from_ice_edge_daily_calipso.py ' +
        '<year> <hemisphere> <dir>' )
    exit()

year = sys.argv[1]
hemisphere = sys.argv[2] 
outpath = sys.argv[3]

if hemisphere[0:2] == 'nh':
    lat0 = 50.
    lat1 = 90.
else:
    lat0 = -90.
    lat1 = -50.    

#read calipso  daily data
path = "/work/Youtong.Zheng/Data/satellite/"
print(path + 'calipso/3D_CFrac/*_' + year + '*.nc')
f = xr.open_mfdataset(glob.glob(path + 'calipso/3D_CFrac/*_' + year + '*.nc'),
                      combine = 'by_coords',decode_cf = True)
f0 = xr.open_mfdataset(glob.glob(path + 'calipso/2D_phase/*_' + year + '*.nc'),
                      combine = 'by_coords',decode_cf = True)

#find indexes of the overlapping time
ind = []
for i, d in enumerate(f.time.values):
    tmp = np.where(f0.time.values == d)[0]
    ind.append(tmp[0])

ind = np.array(ind)

#assign the some variables to the calipso f datsest 
f0 = f0.isel(time = ind)
f0['time'] = f.time
f = f.assign(cllcalipso_ice=f0.cllcalipso_ice).assign(cllcalipso_liq=f0.cllcalipso_liq)

f0 = xr.open_mfdataset(glob.glob(path + 'calipso/*' + year + '*.nc'),
                      combine = 'by_coords',decode_cf = True)

#find indexes of the overlapping time
ind = []
for i, d in enumerate(f.time.values):
    tmp = np.where(f0.time.values == d)[0]
    ind.append(tmp[0])

ind = np.array(ind)

#assign the some variables to the calipso f datsest 
f0 = f0.isel(time = ind)
f0['time'] = f.time
f = f.assign(cllcalipso=f0.cllcalipso)

#some preparing works
f = f.rename({'latitude':'lat'}).rename({'longitude':'lon'})
f = f.sel(lat = slice(lat0, lat1))
lon_plot = f.lon.values
lat_plot = f.lat.values

#read sea ice daily data
f0 = xr.open_mfdataset(glob.glob(path + 'seaice/*' + hemisphere[0:2] + '*' + year + '*.nc'),
                      combine = 'by_coords',decode_cf = True)
f0 = f0.rename({'tdim': 'time'})
f0 = f0.interp(lat=lat_plot, lon=lon_plot)

#find indexes of the overlapping time
ind = []
for i, d in enumerate(f.time.values):
    tmp = np.where(f0.time.values == d)[0]
    ind.append(tmp[0])

ind = np.array(ind)

#print out to double-check
print(f.time.values)
print(ind)

#assign the seaice to the calipso f datsest 
f0 = f0.isel(time = ind).cdr_seaice_conc_gridded
f0['time'] = f.time
f = f.assign(seaice_conc=f0)

#read ERA5 daily
print(glob.glob('/archive/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr-timestep/annual_file-range/Wind/v_10m/*' + year + '.nc'))
f0 = xr.open_dataset(glob.glob('/archive/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr-timestep/annual_file-range/Wind/v_10m/*' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0 = f0.rename({'latitude':'lat'}).rename({'longitude':'lon'})
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).v10
f0['time'] = f.time
f = f.assign(v10=f0)

f0 = xr.open_dataset(glob.glob('/archive/uda/ERA5/Hourly_Data_On_Single_Levels/reanalysis/global/1hr-timestep/annual_file-range/Wind/u_10m/*' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0 = f0.rename({'latitude':'lat'}).rename({'longitude':'lon'})
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).u10
f0['time'] = f.time
f = f.assign(u10=f0)

#Read NCEP profile
f0 = xr.open_dataset(glob.glob('/archive/uda/NCEP-Reanalysis/air.' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).air
f0['time'] = f.time
f = f.assign(air=f0)

f0 = xr.open_dataset(glob.glob('/archive/uda/NCEP-Reanalysis/shum.' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).shum
f0['time'] = f.time
f = f.assign(shum=f0)

f0 = xr.open_dataset(glob.glob('/archive/uda/NCEP-Reanalysis/rhum.' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).rhum
f0['time'] = f.time
f = f.assign(rhum=f0)

f0 = xr.open_dataset(glob.glob('/archive/uda/NCEP-Reanalysis/omega.' + year + '.nc')[0])
f0 = f0.resample(time='1D').mean()
f0.coords['lon'] = (f0.coords['lon'] + 180) % 360 - 180
f0 = f0.sortby(f0.lon).sortby(f0.lat)
f0 = f0.interp(lat=lat_plot, lon=lon_plot)
f0 = f0.isel(time = ind).omega
f0['time'] = f.time
f = f.assign(omega=f0)

#Extract values for computation
# nh and sh corresponds to the default regions of arctic atlantic and antarctic weldell sea shown in the paper
# nh-chukchi and sh-ross corresponds to the arctic chukchi sea and antarctic ross sea, respectively

if hemisphere == 'nh':
    f0 = f.sel(lat=slice(lat0, lat1)).sel(lon=slice(-15, 55))
elif hemisphere == 'nh-chukchi':
    f0 = xr.concat([f.sel(lat=slice(lat0, lat1)).sel(lon=slice(-180, -165)),
                    f.sel(lat=slice(lat0, lat1)).sel(lon=slice(170, 180))], dim="lon")
elif hemisphere == 'sh':
    f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-50, 0))
else:
    f0 = f.sel(lat=slice(-85, -50)).sel(lon=slice(-180, -130))

v_achv = f0.v10.values
u_achv = f0.u10.values
clice_achv = f0.cllcalipso_ice.values
clliq_achv = f0.cllcalipso_liq.values
cll_achv = f0.cllcalipso.values
ice_achv = f0.seaice_conc.values
cl_achv = f0.clcalipso.values
air_achv = f0.air.values
shum_achv = f0.shum.values
rhum_achv = f0.rhum.values
omega_achv = f0.omega.values

lat = f0.lat.values
lon = f0.lon.values
time = f0.time.values 
print(f0.isel(time = 0).alt_mid.values)
altitude = f0.isel(time = 0).alt_mid.values
if hemisphere == 'nh-chukchi':
    altitude = f0.isel(time = 0).alt_mid.values[0]
level = f0.level.values
nlat = len(lat)
nlon = len(lon)
ntime = len(time)
naltitude = len(altitude)
nlevel = len(level)

#computing
print("Start extracting")
print(cl_achv.shape)

cl_out, ice_out = Binning_distance_from_ice_edge(ice_achv, cl_achv, ntime, nlat, nlon, hemisphere[0:2])
air_out, ice_out = Binning_distance_from_ice_edge(ice_achv, air_achv, ntime, nlat, nlon, hemisphere[0:2])
shum_out, ice_out = Binning_distance_from_ice_edge(ice_achv, shum_achv, ntime, nlat, nlon, hemisphere[0:2])
rhum_out, ice_out = Binning_distance_from_ice_edge(ice_achv, rhum_achv, ntime, nlat, nlon, hemisphere[0:2])
omega_out, ice_out = Binning_distance_from_ice_edge(ice_achv, omega_achv, ntime, nlat, nlon, hemisphere[0:2])
v_out, ice_out, count_out = Binning_distance_from_ice_edge(ice_achv, v_achv, ntime, nlat, nlon, hemisphere[0:2], outputcount = True)
u_out, ice_out, count_out = Binning_distance_from_ice_edge(ice_achv, u_achv, ntime, nlat, nlon, hemisphere[0:2], outputcount = True)
clice_out, ice_out, count_out = Binning_distance_from_ice_edge(ice_achv, clice_achv, ntime, nlat, nlon, hemisphere[0:2], outputcount = True)
clliq_out, ice_out, count_out = Binning_distance_from_ice_edge(ice_achv, clliq_achv, ntime, nlat, nlon, hemisphere[0:2], outputcount = True)
cll_out, ice_out, count_out = Binning_distance_from_ice_edge(ice_achv, cll_achv, ntime, nlat, nlon, hemisphere[0:2], outputcount = True)

#make the ds to output
ds_tmp = xr.DataArray(cl_out, coords=[altitude, np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["altitude", "distance","time"]).rename('clcalipso')

ds_tmp1 = xr.DataArray(ice_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('seaice_conc')

ds_tmp2 = xr.DataArray(v_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('v10')

ds_tmp2_add = xr.DataArray(u_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('u10')

ds_tmp3 = xr.DataArray(clice_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('cllcalipsoice')

ds_tmp4 = xr.DataArray(clliq_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('cllcalipsoliq')

ds_tmp5 = xr.DataArray(count_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('count')

ds_tmp6 = xr.DataArray(air_out, coords=[level,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level","distance","time"]).rename('air')

ds_tmp7 = xr.DataArray(shum_out, coords=[level,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level","distance","time"]).rename('shum')

ds_tmp8 = xr.DataArray(rhum_out, coords=[level,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level","distance","time"]).rename('rhum')

ds_tmp9 = xr.DataArray(omega_out, coords=[level,np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["level","distance","time"]).rename('omega')

ds_tmp10 = xr.DataArray(cll_out, coords=[np.arange(-nlat + 0.5, nlat + 0.5, 1.), f0.time.values],
                      dims=["distance","time"]).rename('cllcalipso')

ds_out = xr.merge([ds_tmp, ds_tmp1, ds_tmp2, ds_tmp2_add, ds_tmp3, ds_tmp4, ds_tmp5, ds_tmp6, ds_tmp7, ds_tmp8, ds_tmp9, ds_tmp10]) 

ds_out.attrs['title'] = 'year:' + year

# outpath = '/nbhome/Youtong.Zheng/'
ds_out.to_netcdf(outpath + hemisphere + '_' + year + "_daily_transect_calipso.nc")
print('Job accomplished for year:' + year)
