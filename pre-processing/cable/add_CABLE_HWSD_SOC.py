import numpy as np
import netCDF4 as nc
import xarray as xr



ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables.nc", decode_times=False)

SOC_data = xr.DataArray(
    np.full((7,112,192), -99999.0),
    dims=('soil_hwsd','lat', 'lon'),
    coords={'soil_hwsd': np.arange(7),'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': "CABLE HWSD SOC",
        'units': 'g/kg',
    }
)

HWSD_BD = xr.DataArray(
    np.full((7,112,192), -99999.0),
    dims=('soil_hwsd','lat', 'lon'),
    coords={'soil_hwsd': np.arange(7),'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': "HWSD bulk density",
        'units': 'kg/m3',
    }
)

for i in range(1,8):
    print(i)
    dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_content_new_cable_D"+str(i)+".nc", decode_times=False)
    SOC=dataset['Band1'][:]
    SOC = SOC.fillna(-99999)
    SOC_data[i-1,:]=SOC.values
    
    ds_bulk_density = xr.open_dataset("/data1/zxy/SOC_data_calibration/HWSD/HWSD_SOC_bulk_cable_D"+str(i)+".nc")
    bulk_density=ds_bulk_density['Band1']*1000
    bulk_density = bulk_density.fillna(-99999)
    HWSD_BD[i-1,:]=bulk_density.values    
    

ds['HWSD_SOC']=SOC_data
ds['HWSD_bulk_density']=HWSD_BD

ds.to_netcdf("/data1/zxy/SOC_data_calibration/cable_data/result/CABLE_Time_invariant_variables_new.nc")
