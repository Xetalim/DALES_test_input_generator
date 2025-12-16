import warnings
import pathlib
import logging
from helper_scripts.logging_wrapper import logwrap


def require_condition(expression, error_msg, warn_instead=False):
    if not expression:
        if not warn_instead:
            raise ValueError(error_msg)
        else:
            warnings.warn(error_msg)


def check_ibm_input_settings(config, machine_conf, nml):
    example = "See cases/IBM_small_rect_building.yaml as an example"
    if ("ibm_modifications" in config) or ("NAMIBM" in nml):
        require_condition(
            "NAMIBM" in nml,
            f"Must define NAMIBM to use ibm_modifications. {example}",
        )
        require_condition(
            "ibm_modifications" in config,
            f"Must define ibm_modifications to use NAMIBM. {example}",
        )
    else:
        return
    # assume we are using IBM, check the iadv values
    if ("ibm_modifications" in config) and ("NAMIBM" in nml):
        require_condition(
            nml["NAMIBM"]["lapply_ibm"],
            "lapply_ibm is off but ou request ibm_modifications",
        )


#       NAMIBM:
#     lapply_ibm: true # Activate immersed boundary method
#     lwallheat: false # No heat flux from walls
#     thlibm: 301 # IBM reference theta_l [K] (not used)
#     thlwall: 301 # Wall theta_l [K] (not used)
#     thlroof: 301 # Roof theta_l [K] (not used)
#     qtibm: 0.001 # IBM reference qt [kg/kg]
#   DYNAMICS:
#     iadv_mom: 2 # necessary for IBM to work
#     iadv_tke: 2 # necessary for IBM to work
#     iadv_thl: 2 # necessary for IBM to work
#     iadv_qt: 2 # necessary for IBM to work
#     iadv_sv: 2 # necessary for IBM to work
#     ibas_prf: 2 # necessary for IBM to work
def check_emission_input_settings(
    config, machine_conf, nml, required_files, required_folder_list
):
    if "emissions" in config:
        required_folder_list.append("emissions")


def check_lsm_input_settings(config, machine_conf, nml, required_files):
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


def check_radiation_required_files(machine_conf, nml, exp_id, required_files):
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


def check_config_contains_origin(config):
    if not ("x0" in config["profile"]["grid"]):
        warnings.warn("No origin x0 found in namelist, setting to 0")
        config["profile"]["grid"]["x0"] = 0
    if not ("y0" in config["profile"]["grid"]):
        config["profile"]["grid"]["y0"] = 0
        warnings.warn("No origin y0 found in namelist, setting to 0")


def validate_config_paths(config, machine_conf):
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

    return output_path


def check_core_amounts(machine_conf, nml):
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
