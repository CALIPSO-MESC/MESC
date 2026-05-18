import numpy as np
import netCDF4 as nc
import xarray as xr


dataset=xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/NPP/mean_annual_npp1990-2009.nc", decode_times=False)
ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant.nc", decode_times=False)

# FluxCtolitter(time, mlitter, land)
npp=dataset['npp'][0,:]
npp_filled = npp.fillna(-99999)
print(npp.shape)
npp_new = xr.DataArray(
    npp_filled,
    dims=('lat', 'lon'),
    coords={'lat': ds['lat'], 'lon': ds['lon']},
    attrs={
        'long_name': 'mean annual NPP',
        'units': 'kg C m-2 s-1',
    }
)

ds['npp']=npp_new


ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new.nc")

