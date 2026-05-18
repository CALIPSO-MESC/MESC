import netCDF4 as nc
import pandas as pd
import numpy as np
import xarray as xr

ds = xr.open_dataset("/data1/zxy/SOC_data_calibration/ORCHIDEE_data/data_extraction_ORCHIDEE_model_calibration.nc", decode_times=False)
print(ds.max())
print(ds)

for var_name in ds.data_vars:
    unit = ds[var_name].attrs.get('units', ds[var_name].attrs.get('unit', '-'))
    print(f"Variable: {var_name}({unit})")
    print(f"  Min: {ds[var_name].min().values}")
    print(f"  Max: {ds[var_name].max().values}")
    print("-" * 30)  # 分隔线











