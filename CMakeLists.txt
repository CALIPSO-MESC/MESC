cmake_minimum_required(VERSION 3.18)
project(MESC LANGUAGES Fortran)

# ------------------------------------------------------------
# 0) Fortran standard
# ------------------------------------------------------------
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# 1) Let CMake find our custom modules
# ------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ------------------------------------------------------------
# 2) Find NetCDF Fortran
#   - Prefer nf-config from module environment
#   - Fallback to NetCDFFortran_ROOT / env NETCDF_FORTRAN / env NETCDF
# ------------------------------------------------------------
find_package(NetCDF REQUIRED)
find_package(NetCDFFortran REQUIRED)

message(STATUS "NetCDFFortran include: ${NetCDFFortran_INCLUDE_DIRS}")
message(STATUS "NetCDFFortran libs   : ${NetCDFFortran_LIBRARIES}")

# ------------------------------------------------------------
# 3) Sources
#   - `mesc` static library holds the model modules
#   - `main` is the standalone MESC driver
#   - `mock_lsm_driver` exercises the coupling API in place of an LSM
# ------------------------------------------------------------
add_library(mesc STATIC
  src/mod_precision.f90
  src/mod_namelist.f90
  src/mod_constants.f90
  src/mod_variables.f90
  src/mod_model_core.f90
  src/mod_inout.f90
  src/mod_interface.f90
  src/mod_coupling.f90
  src/mod_calcost.f90
  src/mod_functions.f90
)

# Put .mod files into build/mod
set_target_properties(mesc PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod"
)

# ------------------------------------------------------------
# 4) Include + Link
#   - include: netcdf.mod location
#   - link: nf-config --flibs output or fallback libs
#   (PUBLIC so that executables inherit them)
# ------------------------------------------------------------
target_include_directories(mesc PUBLIC
  "${CMAKE_BINARY_DIR}/mod"
  ${NetCDF_INCLUDES}
  ${NetCDFFortran_INCLUDE_DIRS}
)

# NetCDFFortran_LIBRARIES may be a list of raw flags (e.g., -L... -lnetcdff -lnetcdf ...)
# CMake is OK with passing these through to the linker.
target_link_libraries(mesc PUBLIC
  ${NetCDF_LIBRARIES}
  ${NetCDFFortran_LIBRARIES}
)

# ------------------------------------------------------------
# 5) Intel tweaks (optional)
# ------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(mesc PRIVATE -qopenmp -O3 -fp-model precise)
  target_link_options(mesc PUBLIC -qopenmp)
endif()

# ------------------------------------------------------------
# 6) Executables
#   (name `main` matches your previous -o main)
# ------------------------------------------------------------
add_executable(main src/main.f90)
target_link_libraries(main PRIVATE mesc)

# Mock LSM driver: tests the coupling API end-to-end (runs from test/)
add_executable(mock_lsm_driver test/mock_lsm_driver.f90)
target_link_libraries(mock_lsm_driver PRIVATE mesc)

# ------------------------------------------------------------
# 7) POST-BUILD: copy executables to test directory
# ------------------------------------------------------------
add_custom_command(TARGET main POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
          $<TARGET_FILE:main>
          ${CMAKE_SOURCE_DIR}/test/main
  COMMENT "Copying `main` executable to test directory"
)

add_custom_command(TARGET mock_lsm_driver POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
          $<TARGET_FILE:mock_lsm_driver>
          ${CMAKE_SOURCE_DIR}/test/mock_lsm_driver
  COMMENT "Copying `mock_lsm_driver` executable to test directory"
)
