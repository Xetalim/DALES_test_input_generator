# Dales input case generator
This repository contains a script (main.py) that generates complete DALES cases (except for the dales binary).
Each case is configurable with a yaml config file found in the [cases](./cases) folder.

## Usage
Edit the `machine_conf.yaml` file to point to your DALES source, to your DALES binary.
The `BASE_OUTPUT_PATH` variable in the configuration is the base path of each case, but you can edit the specific path of each case individually in each case yaml file. Tildes (~) do not work in this path (at least on my device).
The required packages can be found in requirements.txt.

You can run the script as a CLI tool where you give the path to the case as an argument to the script.
In each case folder will be a job.001 file which is the job that runs DALES. In the generated input folder will be all input files.


## Changelog:
 -Allowed RTE-RRTMGP runs to be done. Added machine_conf.yaml to prevent having to edit scripts to customise.
 -Added some SLURB model input files (not yet useful)