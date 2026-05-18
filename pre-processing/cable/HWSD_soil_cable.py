import xarray as xr
import numpy as np


ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

bulk_density_raw = ds['HWSD_bulk_density']   # (soil_hwsd, lat, lon)
soc_raw = ds['HWSD_SOC']                     # (soil_hwsd, lat, lon)

bulk_density_new = bulk_density_raw.rename({'soil_hwsd': 'ms'})
soc_new = soc_raw.rename({'soil_hwsd': 'ms'})

bulk_density_new.name = 'Bulk density'
soc_new.name = 'SOC concentration'

soil_layer = xr.DataArray(
    ['0-20cm', '20-40cm', '40-60cm', '60-80cm', '80-100cm', '100-150cm', '150-200cm'],
    dims='ms',
    coords={'ms': np.arange(7)},
)


#=============================================Clay and Silt===========================================#
Clay_fraction = xr.full_like(bulk_density_new, -99999)
Silt_fraction = xr.full_like(bulk_density_new, -99999)
gravel_fraction = xr.full_like(bulk_density_new, -99999)
ph_data       = xr.full_like(bulk_density_new, -99999)

Clay_fraction.name = 'Clay fraction'
Clay_fraction.attrs['units'] = '%'  
Clay_fraction.attrs['long_name'] = 'Clay fraction'

Silt_fraction.name = 'Silt fraction'
Silt_fraction.attrs['units'] = '%'  
Silt_fraction.attrs['long_name'] = 'Silt fraction'

gravel_fraction.name = 'gravel fraction'
gravel_fraction.attrs['units'] = '%'
gravel_fraction.attrs['long_name'] = 'gravel fraction'

ph_data.attrs['long_name'] = 'ph'


for i in range(1,8):
    print(i)
    Clay_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_Clay_cable_D"+str(i)+".nc", decode_times=False)
    Clay=Clay_dataset['Band1'][:]
    Clay = Clay.fillna(-99999)
    Clay_fraction[i-1,:]=Clay.values
    
    Silt_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_Silt_cable_D"+str(i)+".nc", decode_times=False)
    Silt=Silt_dataset['Band1'][:]
    Silt = Silt.fillna(-99999)
    Silt_fraction[i-1,:]=Silt.values
    
    ph_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_ph_cable_D"+str(i)+".nc", decode_times=False)
    ph=ph_dataset['Band1'][:]
    ph = ph.fillna(-99999)
    ph_data[i-1,:]=ph.values    

    gravel_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_coarse_cable_D"+str(i)+".nc", decode_times=False)
    gravel=gravel_dataset['Band1'][:]
    gravel = gravel.fillna(-99999)
    gravel_fraction[i-1,:]=gravel.values



#=================================================Alo==================================================#

Ald_data = xr.DataArray(
    np.full((5,112,192), -99999.0),
    dims=('ms1','lat', 'lon'),
    coords={'ms1': np.arange(5),'lat': ds['lat'], 'lon': ds['lon']},)
Alo_data = xr.DataArray(
    np.full((5,112,192), -99999.0),
    dims=('ms1','lat', 'lon'),
    coords={'ms1': np.arange(5),'lat': ds['lat'], 'lon': ds['lon']},)
Fed_data = xr.DataArray(
    np.full((5,112,192), -99999.0),
    dims=('ms1','lat', 'lon'),
    coords={'ms1': np.arange(5),'lat': ds['lat'], 'lon': ds['lon']},)
Feo_data = xr.DataArray(
    np.full((5,112,192), -99999.0),
    dims=('ms1','lat', 'lon'),
    coords={'ms1': np.arange(5),'lat': ds['lat'], 'lon': ds['lon']},)
    

for i in range(0,5):
    print(i)
    Ald_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ren-gcb2024/reactive-minerals/Ald_"+str(i*20)+"_"+str(i*20+20)+"_cable.nc", decode_times=False)
    Ald=Ald_dataset['Band1'][:]
    Ald = Ald.fillna(-99999)
    Ald_data[i-1,:]=Ald.values
    
    Alo_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ren-gcb2024/reactive-minerals/Alo_"+str(i*20)+"_"+str(i*20+20)+"_cable.nc", decode_times=False)
    Alo=Alo_dataset['Band1'][:]
    Alo = Alo.fillna(-99999)
    Alo_data[i-1,:]=Alo.values
    
    Fed_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ren-gcb2024/reactive-minerals/Fed_"+str(i*20)+"_"+str(i*20+20)+"_cable.nc", decode_times=False)
    Fed=Fed_dataset['Band1'][:]
    Fed = Fed.fillna(-99999)
    Fed_data[i-1,:]=Fed.values
    
    Feo_dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ren-gcb2024/reactive-minerals/Feo_"+str(i*20)+"_"+str(i*20+20)+"_cable.nc", decode_times=False)
    Feo=Feo_dataset['Band1'][:]
    Feo = Feo.fillna(-99999)
    Feo_data[i-1,:]=Feo.values    

soil_layer1 = xr.DataArray(
    ['0-20cm', '20-40cm', '40-60cm', '60-80cm', '80-100cm'],
    dims='ms1',
    coords={'ms1': np.arange(5)},
)

#bulk_density_new.attrs['units'] = 'kg/m3'   
#soc_new.attrs['units'] = 'g/kg'             


ds_new = xr.Dataset({
      'Bulk density': bulk_density_new,
      'SOC concentration': soc_new,
      'Clay fraction': Clay_fraction,
      'Silt fraction': Silt_fraction,
      'gravel fraction': gravel_fraction,
      'ph': ph_data,      
      'Ald': Ald_data,
      'Alo': Alo_data,
      'Fed': Fed_data,
      'Feo': Feo_data,
})

ds_new.attrs['Soil layers (ms=7)'] = ["0-20cm", "20-40cm", "40-60cm", "60-80cm", "80-100cm", "100-150cm", "150-200cm"]
ds_new.attrs['Soil layers1 (ms1=5)'] = ["0-20cm", "20-40cm", "40-60cm", "60-80cm", "80-100cm"]

ds_new.to_netcdf('/data1/zxy/SOC_data_calibration/cable_data/result/hwsd-soil-cable.nc',encoding={
        'Bulk density': {'_FillValue': -99999.0}, 
        'SOC concentration': {'_FillValue': -99999.0}, 
        'Clay fraction': {'_FillValue': -99999.0},
        'Silt fraction': {'_FillValue': -99999.0},
        'gravel fraction': {'_FillValue': -99999.0},
        'ph': {'_FillValue': -99999.0},  
        'Ald': {'_FillValue': -99999.0},   
        'Alo': {'_FillValue': -99999.0},  
        'Fed': {'_FillValue': -99999.0},        
        'Feo': {'_FillValue': -99999.0},
    })








