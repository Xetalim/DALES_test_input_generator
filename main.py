import numpy as np
import os
import glob
import f90nml
import yaml #install pyyaml
import subprocess
import pathlib
import re
# Custom Python scripts/tools/...
from helper_scripts import do_profiles
from helper_scripts import landuse_types
from helper_scripts import create_lsm
from helper_scripts import geometry_modification
def apply_job_conf(job_conf, output_filepath, required_files):
    with open("input_template/job.001", "r") as f:
        content = f.read()
    
    required_transfers = ""
    if len(required_files) != 0:
        for file, originalpath in required_files.items():
            # required_transfers += f'rsync -a "{originalpath}" "$RUNDIR/{file}" || exit\n'
            required_transfers += f'ln -sf "$(pwd)/input/{file}" "$RUNDIR/{file}" || exit\n'
    
    content = content.replace("{{required_filetransfers}}", required_transfers)
    for varname, replacement in job_conf.items():
        content = content.replace(f"{{{{{varname}}}}}", f'"{str(replacement)}"')
    print(f"Writing job file to {output_filepath}")
    with open(output_filepath, "w") as f:
        f.write(content)
def generate_case(config, machine_conf):
    if (config["output"]["path"] is None) and (machine_conf["case_conf"]["BASE_OUTPUT_PATH"] is None):
        raise ValueError(f"Require path in case YAML file or in the machine conf yaml file to output DALES files. (choose a valid path)\nInput files will be created in USER_GIVEN_PATH/{config['output']['name']}/input")
    else:
        if machine_conf["case_conf"]["BASE_OUTPUT_PATH"]:
            output_path = pathlib.Path(machine_conf["case_conf"]["BASE_OUTPUT_PATH"]) / config["output"]["name"]
        else:
            output_path = pathlib.Path(config["output"]["path"]) / config["output"]["name"]
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path / "input", exist_ok=True)

    namoptions = pathlib.Path.cwd() / "input_template" / "input" / "namoptions.001"
    nml = f90nml.read(namoptions)


    # Override namelist parameters from YAML
    for section, values in config["namelist_overrides"].items():
        print(f"Overwriting namelist variables in section {section}")
        if not (section in nml.keys()):
            nml[section] = {}
        for key, value in values.items():
            print(f"\t{key} = {value}")
            nml[section][key] = value
    

    exp_id = nml["RUN"]["iexpnr"]
    required_files = {}
    if nml["PHYSICS"]["iradiation"] == 1:
        required_files[f"backrad.inp.{exp_id:03d}"] = (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "cases" / "example" / "backrad.inp.001").as_posix()
    if nml["PHYSICS"]["iradiation"] == 4:
        # this is an RRTMG case, we need RRTMG_LW and RRTMG_SW and backrad.inp.001.nc
        required_files[f"backrad.inp.{exp_id:03d}.nc"] = (pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc").as_posix()
        required_files[f"rrtmg_lw.nc"] = (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "external" / "RRTMG" / "RRTMG_LW" / "data" / "rrtmg_lw.nc").as_posix()
        required_files[f"rrtmg_sw.nc"] = (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "external" / "RRTMG" / "RRTMG_SW" / "data" / "rrtmg_sw.nc").as_posix()
    elif nml["PHYSICS"]["iradiation"] == 5:
        # this is RTE_RRTMG, we need all data from RTE_RRTMG
        required_files[f"backrad.inp.{exp_id:03d}.nc"] = (pathlib.Path.cwd() / "extra_data/backrad.inp.001.nc").as_posix()
        required_files[f"rrtmg_lw.nc"] = (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "external" / "RRTMG" / "RRTMG_LW" / "data" / "rrtmg_lw.nc").as_posix()
        required_files[f"rrtmg_sw.nc"] = (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "external" / "RRTMG" / "RRTMG_SW" / "data" / "rrtmg_sw.nc").as_posix()
        for file in (pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "external" / "rrtmgp-data").glob("*.nc"):
            required_files[file.name] = file.as_posix()
    
    if nml["NAMSURFACE"]["isurf"] == 11:
        required_files[f"van_genuchten_parameters.nc"] = pathlib.Path(machine_conf["case_conf"]["SOURCE_PATH"]) / "data" / "van_genuchten_parameters.nc"
    



    do_profiles.output_profiles(config["profile"], nml, output_path / "input", bool(config["output"]["plot"]))


    domain = {'x0'     : 0,
              'y0'     : 0,
              'itot'   : nml['DOMAIN']['itot'],
              'jtot'   : nml['DOMAIN']['jtot'],
              'dx'     : nml['DOMAIN']['xsize'] /  nml['DOMAIN']['itot'],
              'dy'     : nml['DOMAIN']['ysize'] /  nml['DOMAIN']['jtot'],
            }
    
    create_lsm.process_input(landuse_types.lu_types_depac,
                              domain,
                              output_path / "input",
                              exp_id,
                              lplot=True,
                              modify_func=lambda lu_types, lu_dict, lsm_input: geometry_modification.lsm_modify_func(config, lu_types, lu_dict, lsm_input))
    
    if "ibm_modifications" in config:
        x, y, itot, jtot = create_lsm.generate_dales_domain(domain)
        ibm_generator = geometry_modification.ibmCreatorClass(x, y)
        for modification in config["ibm_modifications"]:
            ibm_generator.parse_yaml_name(modification)
        ibm_generator.output_nc(output_path / 'input' / f"ibm.inp_{exp_id:03d}.nc")
    

    apply_job_conf(machine_conf["job_conf"], output_path / "job.001", required_files=required_files)
    for file, originalpath in required_files.items():
        subprocess.call(["rsync", "-vt",originalpath , output_path / f"input/{file}"])
    # subprocess.call(["rsync", "-a",pathlib.Path.cwd() / "job.001", output_path])
    subprocess.call(["rsync", "-a",(pathlib.Path.cwd() / "input_template/input").as_posix() + "/" , output_path / f"input"])
    nml.write(output_path / f"input/namoptions.{exp_id:03d}",force=True)
    print(f"Written input files to {output_path / 'input'}")



if __name__ == "__main__":    
    # Load YAML config
    # to_generate = ["cases/medium_extreme.yaml"]
    # to_generate = ["cases/IBM_small_rect_building.yaml"]

    with open("machine_conf.yaml", "r") as file:
        machine_conf = yaml.safe_load(file)

    to_generate = ["cases/small_half_forest_half_grass.yaml"]
    for case in to_generate:
        print(f"Processing case {case}")
        with open(case, "r") as f:
            config = yaml.safe_load(f)
            generate_case(config, machine_conf)



