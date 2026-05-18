import numpy as np
import netCDF4 as nc
import xarray as xr

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)
ds['USDA_SoilSuborder'].attrs['long_name'] = "USDA SoilSuborder"
ds = ds.rename({'Soil_texture': 'USDA_Soil_texture_class'})

Soil_texture_class=np.array(ds['USDA_Soil_texture_class'])

bulk_density=np.array([1.5105,	1.5635,	1.5635,	1.4575,	1.431,	1.5105,	1.6165,	1.5105,	1.5635,	1.643,	1.696,	1.643])*1000

Soil_texture_class_clipped = Soil_texture_class.clip(1, 12)
bulk_density_data = bulk_density[Soil_texture_class_clipped.astype(int) - 1]

invalid_mask = (Soil_texture_class < 1) | (Soil_texture_class > 12)
bulk_density_data = xr.where(invalid_mask, -99999, bulk_density_data)

bulk_density_var = xr.DataArray(
    bulk_density_data,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'bulk_density',
        'units': 'kg/m3',
    }
)

ds['bulk_density']=bulk_density_var


ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new.nc")

