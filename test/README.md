## Run Instructions

After building the project, the `main` executable is automatically copied into
this directory - the test directory.

The test directory has a subdirectory `input/`, which holds files used as input
for the model. After cloning the repository, it will include some files, which are used when
running the test suite. This is a small subset of a global data set to keep the
memory footprint of the repository small and keep the run time of the test suite
short.

### (Optional) Step 1: Prepare input NetCDF file

You may wish to run the test suite over a full global data set rather than the
small subset. If so, you will need to download the NC files from [Google Drive](https://drive.google.com/file/d/1RrnM_XRPshyo8hhxZFwIvXbjnEIGjSmI/view?usp=drive_link).
MODIS_npp_2001-2009_0.5d.nc
ORCHIDEE_ave_1990_2020_time_varing_PFTave.nc 
ORCHIDEE_time_invariant.nc 
hwsd-soil-orchidee-cluster.nc 
hwsd_aoc_halfdeg.nc


Put the downloaded NetCDF file into the `test/input` directory.

### Step 2: Execute the program

Once the input file is properly placed, run the startup script directly from the
command line with
```bash
./run_main.sh
```
