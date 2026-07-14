import matplotlib.pyplot as plt
import numpy as np
import pathlib
import xarray as xr
from matplotlib.patches import Patch
import logging
import os
import xarray as xr
from matplotlib.colors import ListedColormap, BoundaryNorm

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)

# Correction factor for aspect ratio of plots
ASPECT_CORR = 2

# @logwrap


#     return lsm_input, lu_types


def plot_lsm_cover(lsm_netcdf_path, plot_base_path):
    """
    Generate some standard plots of Land Surface Model input data

    Parameters
    ----------
    lsm_input : LSM_input_DALES
        Class containing Dales input parameters for all LU types.
    plotvars : list
        List of variables to plot.

    Returns
    -------
    None.

    """

    data_vars = dict()
    # os.makedirs(pathlib.Path(plot_base_path) , exist_ok=True)

    ds = xr.open_dataset(lsm_netcdf_path, engine="netcdf4")
    cover_names = []
    short_names = []
    cover_vars = []

    for luname, lushort in zip(ds["luname"].values, ds["lushort"].values):
        cover_val = ds[f'cover_{lushort.decode("utf-8")}']
        if np.sum(cover_val) > 0:
            cover_names.append(luname.decode("utf-8"))
            short_names.append(lushort.decode("utf-8"))
            cover_vars.append(f'cover_{lushort.decode("utf-8")}')

    lu_ids = np.arange(len(short_names))

    stacked = np.stack([ds[v].values for v in cover_vars], axis=0)

    # Dominant cover index per cell
    dominant = np.argmax(stacked, axis=0)

    # ------------------------------------------------------------------
    # DISCRETE COLORMAP (17 DISTINCT COLORS)
    # ------------------------------------------------------------------
    if "short_grass" in cover_names:  # IFS MAPPING, use nice colors!

        # Full label → color mapping
        label_to_color = {
            "crops_mixed_farming": "#f1c232",
            "irrigated_crops": "#ffd966",
            "short_grass": "#b6d7a8",
            "tall_grass": "#6aa84f",
            "evergreen_shrubs": "#38761d",
            "deciduous_shrubs": "#6d9e4e",
            "evergreen_needleleaf": "#3d940b",
            "deciduous_needleleaf": "#6fc63d",
            "evergreen_broadleaf": "#006400",
            "deciduous_broadleaf": "#274e13",
            "mixed_forest_wood": "#2f6b2f",
            "interrupted_forest": "#7f9c6e",
            "desert": "#e6c07b",
            "semidesert": "#d9b38c",
            "tundra": "#cfe2f3",
            "bogs_marshes": "#76a5af",
            "inland_water": "#3c78d8",
            "ocean": "#073763",
            "water_land_mixtures": "#6fa8dc",
            "ice_caps_glaciers": "#ffffff",
            "urban": "#666666",
            "road": "#000000",
            "water": "#3c78d8",
        }

        # cover_vars = list of variables actually used (order matters!)
        labels = cover_names

        colors = [label_to_color[l] for l in labels]

        cmap = ListedColormap(colors)
        bounds = np.arange(-0.5, len(labels) + 0.5)
        norm = BoundaryNorm(bounds, cmap.N)
    else:
        ncat = len(cover_vars)  # number of categories

        # Use a qualitative matplotlib colormap and sample it
        base_cmap = plt.get_cmap("tab20")  # good for categorical data
        colors = base_cmap(np.linspace(0, 1, ncat))

        cmap = ListedColormap(colors)
        bounds = np.arange(-0.5, ncat + 0.5)
        norm = BoundaryNorm(bounds, cmap.N)

    # ------------------------------------------------------------------
    # COORDINATES (ASSUMES REGULAR GRID)
    # ------------------------------------------------------------------

    x = ds.x.values
    y = ds.y.values

    # ------------------------------------------------------------------
    # PLOT
    # ------------------------------------------------------------------
    plt.ioff()
    fig, ax = plt.subplots(figsize=(8, 6))

    mesh = ax.pcolormesh(
        x - np.min(x), y - np.min(y), dominant, cmap=cmap, norm=norm, shading="nearest"
    )

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Dominant Cover Type")

    # ------------------------------------------------------------------
    # CATEGORICAL COLORBAR
    # ------------------------------------------------------------------

    cbar = plt.colorbar(mesh, ax=ax, ticks=np.arange(len(cover_names)))
    cbar.ax.set_yticklabels(cover_names)
    cbar.set_label("Cover type")

    plt.tight_layout()

    os.makedirs(plot_base_path, exist_ok=True)
    fig.savefig(plot_base_path / f"cover.png", dpi=300)
    plt.close()
    plt.ion()
    return


if __name__ == "__main__":
    # Example usage
    lsm_netcdf_path = pathlib.Path(
        "/Users/andrevanginkel/Documents/40_Input_and_Runs/42_Dales_Cases/42.01_generated_cases/002_LSM_tests/input/lsm.inp_001.nc"
    )
    plot_base_path = pathlib.Path(
        "/Users/andrevanginkel/Documents/40_Input_and_Runs/42_Dales_Cases/42.01_generated_cases/002_LSM_tests/profiles"
    )
    os.makedirs(plot_base_path, exist_ok=True)
    plot_lsm_cover(lsm_netcdf_path, plot_base_path)
