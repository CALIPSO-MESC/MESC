# FindNetCDF.cmake
#
# Exports:
#   NetCDF_FOUND
#   NetCDF_INCLUDE_DIRS
#   NetCDF_LIBRARIES   (may contain raw flags from nc-config)
#
# Strategy:
#   1) Use nc-config (best on HPC modules)
#   2) Fallback to prefix hints (NetCDF_ROOT / NETCDF)

find_program(NC_CONFIG_EXECUTABLE nc-config)

# Collect prefix hints
set(_NETCDF_HINTS "")
if(DEFINED NetCDF_ROOT)
  list(APPEND _NETCDF_HINTS "${NetCDF_ROOT}")
endif()
if(DEFINED ENV{NETCDF})
  list(APPEND _NETCDF_HINTS "$ENV{NETCDF}")
endif()

# ------------------------------------------------------------
# 1) Preferred: nc-config
# ------------------------------------------------------------
if(NC_CONFIG_EXECUTABLE)
  execute_process(COMMAND ${NC_CONFIG_EXECUTABLE} --includedir
                  OUTPUT_VARIABLE _nc_inc
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${NC_CONFIG_EXECUTABLE} --libs
                  OUTPUT_VARIABLE _nc_libs
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Basic sanity
  set(NetCDF_INCLUDE_DIRS "${_nc_inc}")
  set(NetCDF_LIBRARIES "${_nc_libs}")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NetCDF DEFAULT_MSG
    NetCDF_INCLUDE_DIRS NetCDF_LIBRARIES
  )
  return()
endif()

# ------------------------------------------------------------
# 2) Fallback: manual search by prefix hints
# ------------------------------------------------------------
find_path(NetCDF_INCLUDE_DIRS
  NAMES netcdf.mod netcdf.inc netcdf.h
  HINTS ${_NETCDF_HINTS}
  PATH_SUFFIXES include include/netcdf
)

find_library(NetCDF_C_LIB
  NAMES netcdf
  HINTS ${_NETCDF_HINTS}
  PATH_SUFFIXES lib lib64
)

set(NetCDF_LIBRARIES "")
if(NetCDF_C_LIB)
  list(APPEND NetCDF_LIBRARIES ${NetCDF_C_LIB})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF DEFAULT_MSG
  NetCDF_INCLUDE_DIRS NetCDF_LIBRARIES
)
