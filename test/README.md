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
small subset. If so, you will need to download the files named
`data_extraction_cable_model_calibration_cluster_new.nc` 

`ORCHIDEE_time_invarian.nc `

``hwsd_aoc_halfdeg.nc`` 

``hwsd-soil-orchidee-cluster.nc`` 

``ORCHIDEE_ave_1990_2020_time_varing_PFTave.nc``  from [Google Drive](https://drive.google.com/drive/folders/1Od3GMuTzgwdi7wnynJZ9GZjPH768tq2A?usp=sharing).

Put the downloaded NetCDF files into the `test/input` directory, renaming `data_extraction_cable_model_calibration_cluster_new.nc` to
overwrite `data_extraction_cable_model_calibration_cluster_test.nc`.

**Note** : 

`ORCHIDEE_time_invarian.nc `

``hwsd_aoc_halfdeg.nc``

``hwsd-soil-orchidee-cluster.nc``

``ORCHIDEE_ave_1990_2020_time_varing_PFTave.nc``must exist under `test/input` to run global_orchidee.

### Step 2: Execute the program

Once the input file is properly placed, run the startup script directly from the
command line with

```bash
./run_main.sh
```
or
```
sbatch ./run_global_orchidee.sh
```
