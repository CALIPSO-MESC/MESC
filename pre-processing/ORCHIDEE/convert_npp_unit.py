import numpy as np
import netCDF4 as nc
import xarray as xr

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new.nc", decode_times=False)

npp=ds['npp'][:]
npp_unit = ds['npp'].attrs.get('units', None)
print(npp_unit)

if npp_unit == "kg C m-2 s-1":
    print("npp的单位是"+str(npp_unit))
    ds['npp'] = ds['npp'].where(ds['npp'] != -99999) * 365 * 24 * 3600 * 1000
    ds['npp'] = ds['npp'].fillna(-99999)
    ds['npp'].attrs['unit'] = "g C/m2/year"
    ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/ORCHIDEE_time_invariant_new_new.nc")


