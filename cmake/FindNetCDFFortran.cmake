# NetCDFFortran detection module.
#
# Exports:
#   NETCDF_INCLUDE_DIRS  - include dirs for netcdf.mod, netcdf.h
#   NETCDF_LINK_FLAGS    - link flags for both Fortran and C libraries
#
# Priority:
#   (a)  pkg-config (netcdf-fortran + netcdf)
#   (b)  nf-config  (--includedir / --flibs)
#   (c)  FindNetCDF.cmake (respects NETCDF_FORTRAN, NETCDF, etc.)

set(NETCDF_INCLUDE_DIRS "")
set(NETCDF_LINK_FLAGS "")

# --- (a) pkg-config ---
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_search_module(NF netcdf-fortran)
  pkg_search_module(NC netcdf)
  if(NF_FOUND)
    list(APPEND NETCDF_INCLUDE_DIRS "${NF_INCLUDE_DIRS}")
    list(APPEND NETCDF_LINK_FLAGS ${NF_LDFLAGS})
  endif()
  if(NC_FOUND)
    list(APPEND NETCDF_INCLUDE_DIRS "${NC_INCLUDE_DIRS}")
    list(APPEND NETCDF_LINK_FLAGS ${NC_LDFLAGS})
  endif()
endif()

# --- (b) nf-config ---
if(NOT NETCDF_LINK_FLAGS)
  find_program(NF_CONFIG nf-config)
  if(NF_CONFIG)
    execute_process(COMMAND ${NF_CONFIG} --includedir
                    OUTPUT_VARIABLE _nf_inc
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NF_CONFIG} --flibs
                    OUTPUT_VARIABLE _nf_flibs
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    list(APPEND NETCDF_INCLUDE_DIRS "${_nf_inc}")
    separate_arguments(_nf_flibs UNIX_COMMAND)
    list(APPEND NETCDF_LINK_FLAGS ${_nf_flibs})
  endif()
endif()

# --- (c) FindNetCDF.cmake fallback ---
#    Checks NETCDF_DIR, NETCDF_FORTRAN, NETCDF, CPATH, FPATH, etc.
if(NOT NETCDF_LINK_FLAGS)
  set(NETCDF_F90 "YES")
  find_package(NetCDF REQUIRED)
  set(NETCDF_INCLUDE_DIRS "${NETCDF_INCLUDES}")
  set(NETCDF_LINK_FLAGS "${NETCDF_LIBRARIES}")
endif()

# --- Final reporting ---
list(REMOVE_DUPLICATES NETCDF_INCLUDE_DIRS)

if(NOT NETCDF_LINK_FLAGS)
  message(FATAL_ERROR "NetCDF not found")
endif()

message(STATUS "NetCDF include: ${NETCDF_INCLUDE_DIRS}")
message(STATUS "NetCDF link    : ${NETCDF_LINK_FLAGS}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDFFortran
  DEFAULT_MSG NETCDF_INCLUDE_DIRS NETCDF_LINK_FLAGS)
