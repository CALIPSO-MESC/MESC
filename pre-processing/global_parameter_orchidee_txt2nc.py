#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert Orchidee parameter txt file to NetCDF.

Environment:
    python3.10
    # requires: numpy, netCDF4
    # install:  conda install -c conda-forge numpy netCDF4
    #       or: pip install numpy netCDF4

Usage:
    python global_parameter_orchidee_txt2nc.py

Input:
    data.txt   (line 1 = comment, line 2 = header, line 3+ = data)

Output:
    orchidee_params.nc
"""
import os
import numpy as np
from netCDF4 import Dataset

infile = '../test/input/global_parameter_orchidee.txt'
ncfile = os.path.splitext(infile)[0] + '.nc'

# ---- Read header ----
with open(infile) as f:
    header = None
    n_skip = 0
    for line in f:
        n_skip += 1
        s = line.strip()
        if not s:
            continue                          # skip blank lines
        if s.startswith('!') or s.startswith('#'):
            continue                          # skip comment lines (e.g. "! orchidee parameter ...")
        header = s.split()                    # first non-comment, non-blank line is the header
        break

assert header is not None, f'no header found in {infile}'

idx_name = header[0]                          # 'cluster': index column
varnames = header[1:]                         # 14 parameter names
print(f'header  : {header}')
print(f'varnames: {varnames}')

# ---- Read data (skip all preceding lines: comments + header, i.e. n_skip lines) ----
data = np.loadtxt(infile, skiprows=n_skip)

cluster = data[:, 0].astype(np.int32)
values  = data[:, 1:]

assert values.shape[1] == len(varnames), (values.shape, varnames)
print(f'read {len(cluster)} rows x {values.shape[1]} params')

# ---- Write NetCDF ----
ds = Dataset(ncfile, 'w', format='NETCDF4')

ds.createDimension('orchidee', len(cluster))
ds.createDimension('nparam',   len(varnames))

# Index variable
v = ds.createVariable(idx_name, 'i4', ('orchidee',))
v[:] = cluster
v.long_name = 'orchidee cluster index'
v.units = '-'

# Parameter variables
for i, name in enumerate(varnames):
    v = ds.createVariable(name, 'f8', ('orchidee',))
    v[:] = values[:, i]
    v.long_name = name
    v.units = '-'
    v.missing_value = np.float64(1.e20)

# Parameter name table, for cross-checking on the Fortran side
vn = ds.createVariable('param_names', 'S20', ('nparam',))
for i, name in enumerate(varnames):
    vn[i] = name

ds.close()
print('written:', ncfile)