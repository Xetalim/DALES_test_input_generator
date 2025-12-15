import numpy as np
import os
import glob
import f90nml
import yaml  # install pyyaml
import subprocess
import pathlib
import re
import warnings
import argparse
import logging

# Custom Python scripts/tools/...
from helper_scripts.Atmosphere import do_profiles
from helper_scripts.LSM import landuse_types
from helper_scripts.LSM import create_lsm
from helper_scripts.LBC import openboundary
from helper_scripts.Emission import create_emis
from helper_scripts import geometry_modification
from helper_scripts.check_inputs import (
    validate_config_paths,
    check_core_amounts,
    check_lsm_input_settings,
    check_radiation_required_files,
    check_ibm_input_settings,
    check_emission_input_settings,
    check_config_contains_origin,
)
from helper_scripts.grids import GridDales, get_domain_info_nml

logger = logging.getLogger(__name__)


def apply_job_conf(job_conf, output_filepath, required_files, required_folder_list):
    # sets up the jobfile to link the required files to the run directory
    with open("input_template/job.001", "r") as f:
        content = f.read()

    required_transfers = ""
    if len(required_files) != 0:
        for file, originalpath in required_files.items():
            # required_transfers += f'rsync -a "{originalpath}" "$RUNDIR/{file}" || exit\n'
            required_transfers += (
                f'ln -sf "$(pwd)/input/{file}" "$RUNDIR/{file}" || exit\n'
            )

    content = content.replace("{{required_filetransfers}}", required_transfers)

    required_folders = ""
    if len(required_files) != 0:
        for foldername in required_folder_list:
            # required_transfers += f'rsync -a "{originalpath}" "$RUNDIR/{file}" || exit\n'
            required_folders += f'mkdir -p "$RUNDIR/{foldername}" || exit\n'

    content = content.replace("{{required_folders}}", required_folders)

    required_folder_rsyncs = ""

    for foldername in required_folder_list:
        required_folder_rsyncs += f"for inp in input/{foldername}/*\ndo\nrsync -a $inp $RUNDIR/{foldername}/ || exit\ndone\n"
    content = content.replace("{{required_folder_rsyncs}}", required_folder_rsyncs)

    for varname, replacement in job_conf.items():
        content = content.replace(f"{{{{{varname}}}}}", f'"{str(replacement)}"')
    print(f"Writing job file to {output_filepath}")
    with open(output_filepath, "w") as f:
        f.write(content)


def generate_case(config, machine_conf):
    output_path = validate_config_paths(config, machine_conf)

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path / "input", exist_ok=True)

    # default namelist in input_template/input
    namoptions = pathlib.Path.cwd() / "input_template" / "input" / "namoptions.001"
    nml = f90nml.read(namoptions)

    # Overrides namelist parameters from YAML
    override_namelists(config, nml)

    check_config_contains_origin(config)

    check_required_config_fields(config)

    # exp_id is required to give a name to each input file
    exp_id = nml["RUN"]["iexpnr"]

    grid = GridDales(get_domain_info_nml(nml, config))

    # make sure the amount of gridpoints per core is at least 3
    check_core_amounts(machine_conf, nml)

    # sets up files required to run DALES mostly radiation files
    required_files = {}
    required_folder_list = []

    check_radiation_required_files(machine_conf, nml, exp_id, required_files)

    check_lsm_input_settings(config, machine_conf, nml, required_files)

    check_ibm_input_settings(config, machine_conf, nml)

    check_emission_input_settings(
        config, machine_conf, nml, required_files, required_folder_list
    )

    # create atmospheric profiles
    do_profiles.output_profiles(
        config["profile"],
        exp_id,
        grid,
        output_path / "input",
        bool(config["output"]["plot"]),
    )

    do_lsm(grid, config, output_path, nml, exp_id)

    # set up IBM file and output it
    do_ibm(grid, config, output_path, nml, exp_id)

    # set up SLB file and output it
    do_slb(grid, config, output_path, nml, exp_id)

    # set up open boundary condition input files
    do_openbc(grid, config, output_path, nml, exp_id)

    # set up emission data structure
    emissions = create_emis.setup_emissions(grid, config, output_path, nml, exp_id)

    emissions = create_emis.write_emissions_constant(
        emissions, grid, config, output_path, nml, exp_id
    )

    # set up job script and put it in the case directory
    apply_job_conf(
        machine_conf["job_conf"],
        output_path / "job.001",
        required_files=required_files,
        required_folder_list=required_folder_list,
    )

    # put the required files in the case directory
    write_files(output_path, nml, exp_id, required_files)
    print(f"Written input files to {output_path / 'input'}")


def write_files(output_path, nml, exp_id, required_files):
    for file, originalpath in required_files.items():
        subprocess.call(["rsync", "-t", originalpath, output_path / f"input/{file}"])

    subprocess.call(
        [
            "rsync",
            "-a",
            (pathlib.Path.cwd() / "input_template/input").as_posix() + "/",
            output_path / f"input",
        ]
    )

    subprocess.call(
        [
            "chmod",
            "+x",
            output_path / f"job.001",
        ]
    )

    # write the namelist to the case directory
    nml.write(output_path / f"input/namoptions.{exp_id:03d}", force=True)


def override_namelists(config, nml):
    for section, values in config["namelist_overrides"].items():
        if not (section in nml.keys()):
            nml[section] = {}
        for key, value in values.items():
            nml[section][key] = value


def do_slb(grid: GridDales, config, output_path, nml, exp_id):
    if config["generation_settings"]["useslurb"]:
        return
    slb_generator = geometry_modification.slbCreatorClass(grid)
    if "slb_modifications" in config:
        for modification in config["slb_modifications"]:
            slb_generator.parse_yaml_name(modification)
    slb_generator.output_nc(output_path / "input" / f"inslurb.{exp_id:03d}.nc")


def do_openbc(grid: GridDales, config, output_path, nml, exp_id):
    if not config["generation_settings"]["useopenBC"]:
        return
    openboundary.do_openboundary(grid, config, output_path, nml, exp_id)


def do_ibm(grid: GridDales, config, output_path, nml, exp_id):
    if not ("ibm_modifications" in config):
        return
    ibm_generator = geometry_modification.ibmCreatorClass(grid)
    for modification in config["ibm_modifications"]:
        ibm_generator.parse_yaml_name(modification)
    ibm_generator.output_nc(output_path / "input" / f"ibm.inp_{exp_id:03d}.nc")


def do_lsm(grid: GridDales, config, output_path, nml, exp_id):
    if not config["generation_settings"]["uselsm"]:
        return
    lu_types = landuse_types.lu_types_depac
    if "useslurb" in config["generation_settings"]:
        if config["generation_settings"]["useslurb"]:
            lu_types["slb"] = landuse_types.slb
    create_lsm.process_input(
        lu_types,
        grid,
        output_path / "input",
        exp_id,
        lplot=True,
        modify_func=lambda lu_types, lu_dict, lsm_input: geometry_modification.lsm_modify_func(
            config, lu_types, lu_dict, lsm_input, grid
        ),
    )


def check_required_config_fields(config):
    if not ("generation_settings" in config):
        config["generation_settings"] = {}

    for setting in ["uselsm", "useslurb", "useopenBC"]:
        if not (setting in config["generation_settings"]):
            config["generation_settings"][setting] = False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create input files from yaml file")
    parser.add_argument(
        "casefile",
        help="Path of yaml file (e.g. 'cases/mini.yaml')",
    )
    logging.basicConfig(level=logging.INFO)
    args = parser.parse_args()
    casefile = args.casefile

    # casefile = "my_own_cases/mini_open.yaml"
    # casefile = "my_own_cases/mini_in_mini_open.yaml"
    logger.info("Processing casefile %s", casefile)
    with open("machine_conf.yaml", "r") as file:
        machine_conf = yaml.safe_load(file)
    with open(casefile, "r") as f:
        config = yaml.safe_load(f)
    generate_case(config, machine_conf)
