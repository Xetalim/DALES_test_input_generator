import os
import f90nml
import yaml  # install pyyaml to make this import work
import subprocess
import pathlib
import argparse
import logging
from dask.distributed import Client
import dask

from helper_scripts.Atmosphere import do_profiles
from helper_scripts.LSM import landuse_types
from helper_scripts.LSM import LSM_output_dales
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


def setup_logging(config_path="logging.yaml"):
    import logging.config

    path = pathlib.Path(config_path)
    if path.exists():
        with path.open() as f:
            user_cfg = yaml.safe_load(f)

        logging.config.dictConfig(user_cfg)
    else:
        logging.basicConfig()


def apply_job_conf(
    job_conf, input_template, output_filepath, required_files, required_folder_list
):
    # sets up the jobfile to link the required files to the run directory
    if "job_file" in input_template:
        job_filename = pathlib.Path("input_template") / input_template["job_file"]
    else:
        job_filename = pathlib.Path("input_template") / "job.001"
    with open(job_filename, "r") as f:
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
        # required_folder_rsyncs += f"for inp in input/{foldername}/*\ndo\nrsync -a $inp $RUNDIR/{foldername}/ || exit\ndone\n"
        required_folder_rsyncs += f"for inp in input/{foldername}/*\ndo\nln -sf $inp '$RUNDIR/{foldername}/$(basename $inp)' || exit|| exit\ndone\n"
    content = content.replace("{{required_folder_rsyncs}}", required_folder_rsyncs)

    for varname, replacement in job_conf.items():
        content = content.replace(f"{{{{{varname}}}}}", f'"{str(replacement)}"')
    logger.info(f"Writing job file to {output_filepath}")
    with open(output_filepath, "w") as f:
        f.write(content)


class DalesCase:
    def __init__(self, config, machine_conf):
        self.config = config
        self.machine_conf = machine_conf

        self.output_path = validate_config_paths(config, self.machine_conf)

        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.output_path / "input", exist_ok=True)

        # default namelist in input_template/input
        self.nml = f90nml.read(
            pathlib.Path.cwd() / "input_template" / "input" / "namoptions.001"
        )

        # Overrides namelist parameters from YAML
        override_namelists(self.config, self.nml)

        check_config_contains_origin(self.config)

        check_required_config_fields(self.config)

        # exp_id is required to give a name to each input file
        self.exp_id = self.nml["RUN"]["iexpnr"]

        self.grid = GridDales(get_domain_info_nml(self.nml, self.config))

        # make sure the amount of gridpoints per core is at least 3
        check_core_amounts(self.machine_conf, self.nml)

        # sets up files required to run DALES; mostly radiation files
        self.required_files = {}
        self.required_folder_list = []

        check_radiation_required_files(
            self.machine_conf, self.nml, self.exp_id, self.required_files
        )

        check_lsm_input_settings(
            self.config, self.machine_conf, self.nml, self.required_files
        )

        check_ibm_input_settings(self.config, self.machine_conf, self.nml)

        check_emission_input_settings(
            self.config,
            machine_conf,
            self.nml,
            self.required_files,
            self.required_folder_list,
        )

    def setup_writers(self):

        # create atmospheric profiles
        self.profilewriter = do_profiles.AtmosphereWriter(self.grid)

        if self.config["generation_settings"]["uselsm"]:
            lu_types = landuse_types.lu_types_depac
            if "useslurb" in self.config["generation_settings"]:
                if self.config["generation_settings"]["useslurb"]:
                    lu_types["slb"] = landuse_types.slb

            self.lsm_writer = LSM_output_dales.LSM_output_dales(
                self.grid,
                lu_types,
            )
        if self.config["generation_settings"]["useibm"]:
            self.ibm_generator = geometry_modification.ibmCreatorClass(self.grid)

        if self.config["generation_settings"]["useslurb"]:
            self.slb_generator = geometry_modification.slbCreatorClass(self.grid)

        # set up emission data structure
        self.emissions_creator = create_emis.setup_emissions(
            self.grid, self.config, self.output_path, self.nml, self.exp_id
        )

        if self.config["generation_settings"]["useopenBC"]:
            self.openbc_creator = openboundary.do_openboundary(self.grid)

    def default_from_yaml(self):

        self.profilewriter.apply_profiles(self.config["profile"])

        if self.config["generation_settings"]["uselsm"]:
            self.lsm_writer.standard_fill_geometry_modification(
                modify_func=lambda lsm_input: geometry_modification.lsm_modify_func(
                    self.config, lsm_input, self.grid
                ),
            )
        if self.config["generation_settings"]["useibm"]:
            for modification in self.config["ibm_modifications"]:
                self.ibm_generator.parse_yaml_name(modification)

        if self.config["generation_settings"]["useslurb"]:
            if "slb_modifications" in self.config:
                for modification in self.config["slb_modifications"]:
                    self.slb_generator.parse_yaml_name(modification)

        logger.warning("Emission is not setting up job properly!!")
        create_emis.write_emissions_constant(
            self.emissions_creator,
            self.grid,
            self.config,
            self.output_path,
            self.nml,
            self.exp_id,
        )
        if self.config["generation_settings"]["useopenBC"]:
            self.openbc_creator.setup(
                self.config, self.output_path, self.nml, self.exp_id
            )

    def prepare_calculations(self):
        if self.config["generation_settings"]["useopenBC"]:
            self.openbc_creator.prepare_calculation(
                self.config, self.output_path, self.nml, self.exp_id
            )

    def write_files(self):
        self.profilewriter.write_profiles(
            output_path=self.output_path / "input", exp_id=self.exp_id
        )
        if self.config["output"]["plot"]:
            self.profilewriter.plot_profiles(
                output_path=self.output_path / "input", exp_id=self.exp_id
            )

        if self.config["generation_settings"]["uselsm"]:
            self.lsm_writer.save_netcdf(self.output_path / "input", self.exp_id)

        if self.config["generation_settings"]["useibm"]:
            self.ibm_generator.output_nc(
                self.output_path / "input" / f"ibm.inp_{self.exp_id:03d}.nc"
            )

        if self.config["generation_settings"]["useslurb"]:
            self.slb_generator.output_nc(
                self.output_path / "input" / f"inslurb.{self.exp_id:03d}.nc"
            )
        if self.config["generation_settings"]["useopenBC"]:
            self.openbc_creator.write_openbcs(self.output_path)

        # set up job script and put it in the case directory
        apply_job_conf(
            self.machine_conf["job_conf"],
            self.config["input_template"],
            self.output_path / "job.001",
            required_files=self.required_files,
            required_folder_list=self.required_folder_list,
        )

        # put the required files in the case directory
        write_files(self.output_path, self.nml, self.exp_id, self.required_files)
        logger.info(f"Written input files to {self.output_path / 'input'}")


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


# def do_openbc(grid: GridDales, config, output_path, nml, exp_id):


def do_ibm(grid: GridDales, config, output_path, nml, exp_id):
    if not ("ibm_modifications" in config):
        return
    ibm_generator = geometry_modification.ibmCreatorClass(grid)
    for modification in config["ibm_modifications"]:
        ibm_generator.parse_yaml_name(modification)
    ibm_generator.output_nc(output_path / "input" / f"ibm.inp_{exp_id:03d}.nc")


def check_required_config_fields(config):
    if not ("generation_settings" in config):
        config["generation_settings"] = {}

    for setting in ["uselsm", "useslurb", "useopenBC", "useibm"]:
        if not (setting in config["generation_settings"]):
            config["generation_settings"][setting] = False


### FOR DEBUGGING, SET THE DASK SCHEDULER TO SINGLE-THREADED
### ELSE LEAVE COMMENTED OUT
# dask.config.set(
#     scheduler="single-threaded"
# )  # overwrite default with single-threaded scheduler


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create input files from yaml file")
    parser.add_argument(
        "casefile",
        help="Path of yaml file (e.g. 'cases/mini.yaml')",
    )
    # client = Client()  # start distributed scheduler locally.

    setup_logging("logging.yaml")

    args = parser.parse_args()
    casefile = args.casefile
    # casefile = "my_own_cases/mini_open.yaml"
    logger.info("Processing casefile %s", casefile)

    with open("machine_conf.yaml", "r") as file:
        machine_conf = yaml.safe_load(file)

    with open(casefile, "r") as f:
        config = yaml.safe_load(f)
        case = DalesCase(config, machine_conf=machine_conf)
        case.setup_writers()
        case.default_from_yaml()
        case.prepare_calculations()
        case.write_files()
