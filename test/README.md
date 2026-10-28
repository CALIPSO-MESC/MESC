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

- `data_extraction_cable_model_calibration_cluster_new.nc`

- `ORCHIDEE_time_invariant.nc`

- `hwsd_aoc_halfdeg.nc`

- `hwsd-soil-orchidee-cluster.nc`

- `ORCHIDEE_ave_1990_2020_time_varing_PFTave.nc` 

from [Google Drive](https://drive.google.com/drive/folders/1Od3GMuTzgwdi7wnynJZ9GZjPH768tq2A?usp=sharing).

Put the downloaded NetCDF files into the `test/input` directory, renaming `data_extraction_cable_model_calibration_cluster_new.nc` to
overwrite `data_extraction_cable_model_calibration_cluster_test.nc`.

**Note** :

- `ORCHIDEE_time_invariant.nc`
- `hwsd_aoc_halfdeg.nc`
- `hwsd-soil-orchidee-cluster.nc`
- `ORCHIDEE_ave_1990_2020_time_varing_PFTave.nc` 

must exist under `test/input` to run global_orchidee.

### Step 2: Execute the program

Once the input file is properly placed, run the startup script directly from the command line:

```bash
./run_main.sh
```

Or submit the job to the Slurm scheduler on an HPC system:

```Shell
sbatch ./run_global_orchidee.slm
```

**Important**
This Slurm script is provided as a template only and will not run unmodified on other HPC systems.
Users **must update cluster-specific configurations** before submission:

1. Slurm job parameters
   - `-p`: partition / queue name
   - `-A`: project or allocation account code
   - Job name, number of tasks and CPU cores according to your resource needs
2. Environment modules
   - The example `module load` commands for the Intel OneAPI compiler and NetCDF library are site-specific.
   - Module names, versions and availability differ across HPC platforms.

Guidance for these configurable options is included as comments at the top of `run_global_orchidee.slm`.
