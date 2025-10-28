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

# Custom Python scripts/tools/...
from helper_scripts import do_profiles
from helper_scripts import landuse_types
from helper_scripts import create_lsm
from helper_scripts import geometry_modification


def set_namelist_override(config, namelist, key, value):
    if not (namelist in config["namelist_overrides"]):
        config["namelist_overrides"][namelist] = {}
    config["namelist_overrides"][namelist][key] = value


def apply_job_conf(job_conf, output_filepath, required_files):
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
    for varname, replacement in job_conf.items():
        content = content.replace(f"{{{{{varname}}}}}", f'"{str(replacement)}"')
    print(f"Writing job file to {output_filepath}")
    with open(output_filepath, "w") as f:
        f.write(content)


def generate_case(config, machine_conf):
    if (config["output"]["path"] is None) and (
        machine_conf["case_conf"]["BASE_OUTPUT_PATH"] is None
    ):
        raise ValueError(
            f"Require path in case YAML file or in the machine conf yaml file to output DALES files. (choose a valid path)\nInput files will be created in USER_GIVEN_PATH/{config['output']['name']}/input"
        )
    else:
        if machine_conf["case_conf"]["BASE_OUTPUT_PATH"]:
            output_path = (
                pathlib.Path(machine_conf["case_conf"]["BASE_OUTPUT_PATH"])
                / config["output"]["name"]
            )
        else:
            output_path = (
                pathlib.Path(config["output"]["path"]) / config["output"]["name"]
            )
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path / "input", exist_ok=True)

    # default namelist in input_template/input
    namoptions = pathlib.Path.cwd() / "input_template" / "input" / "namoptions.001"
    nml = f90nml.read(namoptions)

    # Overrides namelist parameters from YAML
    for section, values in config["namelist_overrides"].items():
        print(f"Overwriting namelist variables in section {section}")
        if not (section in nml.keys()):
            nml[section] = {}
        for key, value in values.items():
            print(f"\t{key} = {value}")
            nml[section][key] = value

    # exp_id is required to give a name to each input file
    exp_id = nml["RUN"]["iexpnr"]

    # make sure the amount of gridpoints per core is at least 3
    ncores = machine_conf["job_conf"]["numcores"]
    itot = nml["DOMAIN"]["itot"]
    jtot = nml["DOMAIN"]["jtot"]
    if not (
        ("nprocx" in nml["RUN"])
        and ("nprocy" in nml["RUN"])
        and ("nprocs" in nml["RUN"])
    ):
        if jtot // (ncores) < 3:
            warnings.warn(
                f"Too many cores to divide {jtot=} by {ncores=}. Are you sure you want to simulate with MPI? Setting ncores to 1"
            )
            ncores = 1
            machine_conf["job_conf"]["numcores"] = ncores

    # sets up files required to run DALES mostly radiation files
    required_files = {}
    if nml["PHYSICS"]["iradiation"] == 1:
        required_files[f"backrad.inp.{exp_id:03d}"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "cases"
            / "example"
            / "backrad.inp.001"
        ).as_posix()
    if nml["PHYSICS"]["iradiation"] == 4:
        # this is an RRTMG case, we need RRTMG_LW and RRTMG_SW and backrad.inp.001.nc
        required_files[f"backrad.inp.{exp_id:03d}.nc"] = (
            pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc"
        ).as_posix()
        required_files[f"rrtmg_lw.nc"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "external"
            / "RRTMG"
            / "RRTMG_LW"
            / "data"
            / "rrtmg_lw.nc"
        ).as_posix()
        required_files[f"rrtmg_sw.nc"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "external"
            / "RRTMG"
            / "RRTMG_SW"
            / "data"
            / "rrtmg_sw.nc"
        ).as_posix()
    elif nml["PHYSICS"]["iradiation"] == 5:
        # this is RTE_RRTMG, we need all data from RTE_RRTMG
        required_files[f"backrad.inp.{exp_id:03d}.nc"] = (
            pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc"
        ).as_posix()
        required_files[f"rrtmg_lw.nc"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "external"
            / "RRTMG"
            / "RRTMG_LW"
            / "data"
            / "rrtmg_lw.nc"
        ).as_posix()
        required_files[f"rrtmg_sw.nc"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "external"
            / "RRTMG"
            / "RRTMG_SW"
            / "data"
            / "rrtmg_sw.nc"
        ).as_posix()
        for file in (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "external"
            / "rrtmgp-data"
        ).glob("*.nc"):
            required_files[file.name] = file.as_posix()

    if nml["NAMSURFACE"]["isurf"] != 11:
        config["generation_settings"]["uselsm"] = False
    if not config["generation_settings"]["uselsm"]:
        # skip LSM generation, make sure isurf != 11 (lsm)
        if nml["NAMSURFACE"]["isurf"] == 11:
            warnings.warn("Can't have isurf==1 AND uselsm=false, setting isurf to 2")
            nml["NAMSURFACE"]["isurf"] = 2

    else:
        nml["NAMSURFACE"]["isurf"] = 11

    # LSM requires van_genuchten soil paramaeters, so we add it to the required files list
    if config["generation_settings"]["uselsm"]:
        required_files[f"van_genuchten_parameters.nc"] = (
            pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"])
            / "data"
            / "van_genuchten_parameters.nc"
        )

    # create atmospheric profiles
    do_profiles.output_profiles(
        config["profile"], nml, output_path / "input", bool(config["output"]["plot"])
    )

    if config["generation_settings"]["uselsm"]:
        lu_types = landuse_types.lu_types_depac
        if "useslurb" in config["generation_settings"]:
            if config["generation_settings"]["useslurb"]:
                lu_types["slb"] = landuse_types.slb
        create_lsm.process_input(
            lu_types,
            create_lsm.get_domain_info_nml(nml),
            output_path / "input",
            exp_id,
            lplot=True,
            modify_func=lambda lu_types, lu_dict, lsm_input: geometry_modification.lsm_modify_func(
                config, lu_types, lu_dict, lsm_input
            ),
        )

    # set up IBM file and output it
    if "ibm_modifications" in config:
        x, y, itot, jtot = create_lsm.generate_dales_domain(
            create_lsm.get_domain_info_nml(nml)
        )
        ibm_generator = geometry_modification.ibmCreatorClass(x, y)
        for modification in config["ibm_modifications"]:
            ibm_generator.parse_yaml_name(modification)
        ibm_generator.output_nc(output_path / "input" / f"ibm.inp_{exp_id:03d}.nc")

    # set up job script and put it in the case directory
    apply_job_conf(
        machine_conf["job_conf"], output_path / "job.001", required_files=required_files
    )

    # put the required files in the case directory
    for file, originalpath in required_files.items():
        subprocess.call(["rsync", "-vt", originalpath, output_path / f"input/{file}"])

    subprocess.call(
        [
            "rsync",
            "-a",
            (pathlib.Path.cwd() / "input_template/input").as_posix() + "/",
            output_path / f"input",
        ]
    )
    # write the namelist to the case directory
    nml.write(output_path / f"input/namoptions.{exp_id:03d}", force=True)
    print(f"Written input files to {output_path / 'input'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create input files from yaml file")
    parser.add_argument("casefile", help="Path of yaml file (e.g. 'cases/mini.yaml')")
    try:
        args = parser.parse_args()
        with open("machine_conf.yaml", "r") as file:
            machine_conf = yaml.safe_load(file)
        with open(args.casefile, "r") as f:
            config = yaml.safe_load(f)
        generate_case(config, machine_conf)
    except argparse.ArgumentError:
        with open("machine_conf.yaml", "r") as file:
            machine_conf = yaml.safe_load(file)

        to_generate = ["cases/samptendtest.yaml"]
        for case in to_generate:
            print(f"Processing case {case}")
            with open(case, "r") as f:
                config = yaml.safe_load(f)
                generate_case(config, machine_conf)
