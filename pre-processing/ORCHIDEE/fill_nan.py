import xarray as xr


ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/soiltemp_daily/soiltemp1990_2009_new.nc", decode_times=False)


ds = ds.fillna(-99999.0)

#for var in ds.variables: 
#    ds[var].attrs['_FillValue'] = -99999.0  


ds.to_netcdf("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/soiltemp_daily/soiltemp1990_2009_new1.nc")