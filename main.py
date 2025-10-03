import numpy as np
import os
import glob
import f90nml
import yaml #install pyyaml
import subprocess
import pathlib
# Custom Python scripts/tools/...
from helper_scripts import do_profiles
from helper_scripts import landuse_types
from helper_scripts import create_lsm
BASE_OUTPUT_PATH = None
class modifierClass:
    # class to edit land use types (lu_types) for the DALES LSM input. Set type using set_type and
    # give a mask with any of the geometry primitives.
    def __init__(self, lu_types, lu_dict, lsm_input):
        self.lu_types = lu_types
        self.lu_dict = lu_dict
        self.lsm_input = lsm_input

        self.meshx, self.meshy = np.meshgrid(lsm_input.x, lsm_input.y)
    
    def returnVars(self):
        return self.lu_types, self.lu_dict, self.lsm_input

    def set_type(self, mask, lu_type):
        if not (lu_type in self.lu_dict.keys()):
            raise KeyError(f"Incorrect lu_type given {lu_type}, {self.lu_dict.keys()}")
        if not (lu_type in self.lu_types.keys()):
            raise KeyError(f"Incorrect lu_type given {lu_type}, {self.lu_types.keys()}")
        self.lu_types[lu_type]["lu_frac"][mask] = 1
        for other_lu_type in self.lu_types.keys():
            if lu_type != other_lu_type:
                self.lu_types[other_lu_type]["lu_frac"][mask] = 0
        self.lu_types[lu_type]["lu_frac"][mask] = 1
        for other_lu_type in self.lu_types.keys():
            if lu_type != other_lu_type:
                self.lu_types[other_lu_type]["lu_frac"][mask] = 0
    
    def allGeometry(self):
        return np.ones_like(self.meshx, dtype=bool)
    def circleGeometry_realspace(self, x0, y0, size):
        return ((self.meshx - x0)**2 + (self.meshy - y0)**2 <= size**2)
    def rectangleGeometry_realspace(self, minx, maxx, miny, maxy):
        return (self.meshx >= minx) & (self.meshx <= maxx) & (self.meshy >= miny) & (self.meshy <= maxy)
    def rectangleGeometry_idxspace(self, minx, maxx, miny, maxy):
        idxmesh, idymesh = np.meshgrid(np.arange(len(self.lsm_input.x)), np.arange(len(self.lsm_input.y)))
        return (idxmesh >= minx) & (idxmesh <= maxx) & (idymesh >= miny) & (idymesh <= maxy)
    def circleGeometry_idxspace(self, idx0, idy0, size):
        idxmesh, idymesh = np.meshgrid(np.arange(len(self.lsm_input.x)), np.arange(len(self.lsm_input.y)))
        return ((idxmesh - idx0)**2 + (idymesh - idy0)**2 <= size**2)
    def parse_yaml_name(self, modification):
        name = modification["geometry"]
        param_dic = modification["params"]
        dic = {"circle_real": lambda: self.circleGeometry_idxspace(**param_dic)
               ,"all":lambda: self.allGeometry()
               ,"rectangle_real":lambda: self.rectangleGeometry_realspace(**param_dic)
               ,"rectangle_idx":lambda: self.rectangleGeometry_idxspace(**param_dic)
               ,"circle_idx":lambda: self.circleGeometry_idxspace(**param_dic)}
        self.set_type(dic[name](), modification["type"]) 

def generate_case(config):
    if (config["output"]["path"] is None) and (BASE_OUTPUT_PATH is None):
        raise ValueError(f"Require path in YAML file or in the variable BASE_OUTPUT_PATH to output DALES files. (choose a valid path)\nInput files will be created in USER_GIVEN_PATH/{config['output']['name']}/input")
    else:
        if BASE_OUTPUT_PATH:
            output_path = pathlib.Path(BASE_OUTPUT_PATH) / config["output"]["name"]
        else:
            output_path = pathlib.Path(config["output"]["path"]) / config["output"]["name"]
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path / "input", exist_ok=True)

    namoptions = pathlib.Path.cwd() / "input_template" / "namoptions.001"
    nml = f90nml.read(namoptions)

    # Override namelist parameters from YAML
    for section, values in config["namelist_overrides"].items():
        print(f"Overwriting namelist variables in section {section}")
        for key, value in values.items():
            print(f"\t{key} = {value}")
            nml[section][key] = value
    do_profiles.output_profiles(config["profile"], nml, output_path / "input", bool(config["output"]["plot"]))


    # dx = xsize / float(itot)
    # dy = ysize / float(jtot)
    exp_id = nml["RUN"]["iexpnr"]
    domain = {'x0'     : 0,
              'y0'     : 0,
              'itot'   : nml['DOMAIN']['itot'],
              'jtot'   : nml['DOMAIN']['jtot'],
              'dx'     : nml['DOMAIN']['xsize'] /  nml['DOMAIN']['itot'],
              'dy'     : nml['DOMAIN']['ysize'] /  nml['DOMAIN']['jtot'],
            }
    def modify_func(lu_types, lu_dict, lsm_input):
        # function that edits the DALES LSM land use types using the modifier class.
        modifier = modifierClass(lu_types, lu_dict, lsm_input)
        for modification in config["land_use_modifications"]:
            modifier.parse_yaml_name(modification)
        # modifier.set_type(modifier.allGeometry(), "grs")
        # modifier.set_type(modifier.rectangleGeometry_idxspace(0,32,0,16), "fbd")
        # modifier.set_type(modifier.circleGeometry_realspace(np.mean(modifier.meshx),np.mean(modifier.meshy),1600),"fbd")
        # modifier.set_type(modifier.circleGeometry_idxspace(16,16,4),"urb")
        return modifier.returnVars()
    
    create_lsm.process_input(landuse_types.lu_types_depac,
                              domain,
                              output_path / "input",
                              exp_id,
                              lplot=False,
                              modify_func=modify_func)
    
    subprocess.call(["rsync", "-a",pathlib.Path.cwd() / "job.001", output_path])
    subprocess.call(["rsync", "-a",(pathlib.Path.cwd() / "input_template").as_posix() + "/" , output_path / f"input"])
    nml.write(output_path / f"input/namoptions.{exp_id:03d}",force=True)
    print(f"Written input files to {output_path / 'input'}")


if __name__ == "__main__":    
    # Load YAML config
    for case in glob.glob("cases/*.yaml"):
        print(f"Processing case {case}")
        with open(case, "r") as f:
            config = yaml.safe_load(f)
            generate_case(config)



