import numpy as np
import netCDF4
class modifierClass:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.meshx, self.meshy = np.meshgrid(x, y)
        self.idxmesh, self.idymesh = np.meshgrid(np.arange(len(x)), np.arange(len(y)))    
    def allGeometry(self):
        return np.ones_like(self.meshx, dtype=bool)
    def circleGeometry_realspace(self, x0, y0, size):
        return ((self.meshx - x0)**2 + (self.meshy - y0)**2 <= size**2)
    def rectangleGeometry_realspace(self, minx, maxx, miny, maxy):
        return (self.meshx >= minx) & (self.meshx <= maxx) & (self.meshy >= miny) & (self.meshy <= maxy)
    def rectangleGeometry_idxspace(self, minx, maxx, miny, maxy):
        return (self.idxmesh >= minx) & (self.idxmesh <= maxx) & (self.idymesh >= miny) & (self.idymesh <= maxy)
    def circleGeometry_idxspace(self, idx0, idy0, size):
        return ((self.idxmesh - idx0)**2 + (self.idymesh - idy0)**2 <= size**2)
    def parse_yaml_name(self, modification):
        param_dic = modification["params"]
        dic = {"circle_real": lambda: self.circleGeometry_idxspace(**param_dic)
               ,"all":lambda: self.allGeometry()
               ,"rectangle_real":lambda: self.rectangleGeometry_realspace(**param_dic)
               ,"rectangle_idx":lambda: self.rectangleGeometry_idxspace(**param_dic)
               ,"circle_idx":lambda: self.circleGeometry_idxspace(**param_dic)}
        
        name = modification["geometry"]
        geometry_function = dic[name]
        geometry = geometry_function()
        
        self.do_modification(geometry, modification)
    
    def do_modification(self, geometry, modification):
        raise NotImplementedError("Can't call do_modification of superclass, call the subclass instead")

class LsmModifier(modifierClass):
    # class to edit land use types (lu_types) for the DALES LSM input. Set type using set_type and
    # give a mask with any of the geometry primitives.
    def __init__(self, lu_types, lu_dict, lsm_input):
        self.lu_types = lu_types
        self.lu_dict = lu_dict
        self.lsm_input = lsm_input
        super().__init__(x=lsm_input.x, y=lsm_input.y)
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
    def do_modification(self, geometry, modification):
        self.set_type(geometry, modification["type"])


class ibmCreatorClass(modifierClass):
    def __init__(self, x, y):
        super().__init__(x=x, y=y)
        self.bc_height = np.zeros_like(self.meshx)
    def do_modification(self, geometry, modification):
        self.bc_height[geometry] = modification["height"]
    def output_nc(self, filename):
        with netCDF4.Dataset(filename, 'w') as nc:

            nc.createDimension('x', len(self.x))
            nc.createDimension('y', len(self.y))

            var_x = nc.createVariable('x', float, 'x')
            var_y = nc.createVariable('y', float, 'y')

            var_x[:] = self.x[:]
            var_y[:] = self.y[:]

            dims = ['y', 'x']
            bc_height_var  = nc.createVariable("bc_height", float, dims)
            bc_height_var[:,:] = self.bc_height[:,:]


def lsm_modify_func(config, lu_types, lu_dict, lsm_input):
    # function that edits the DALES LSM land use types using the modifier class.
    modifier = LsmModifier(lu_types, lu_dict, lsm_input)
    for modification in config["land_use_modifications"]:
        modifier.parse_yaml_name(modification)
    # modifier.set_type(modifier.allGeometry(), "grs")
    # modifier.set_type(modifier.rectangleGeometry_idxspace(0,32,0,16), "fbd")
    # modifier.set_type(modifier.circleGeometry_realspace(np.mean(modifier.meshx),np.mean(modifier.meshy),1600),"fbd")
    # modifier.set_type(modifier.circleGeometry_idxspace(16,16,4),"urb")
    return modifier.returnVars()