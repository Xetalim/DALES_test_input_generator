import matplotlib.pyplot as plt
import numpy as np
import pathlib
import xarray as xr
import logging
import os
import math
from matplotlib.colors import ListedColormap, BoundaryNorm

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)

# Correction factor for aspect ratio of plots
ASPECT_CORR = 2

# @logwrap


#     return lsm_input, lu_types


def _decode_text(val):
    if isinstance(val, bytes):
        return val.decode("utf-8")
    return str(val)


def _plot_2d_field(ax, x, y, data_2d, title, cmap="viridis"):
    mesh = ax.pcolormesh(
        x - np.min(x),
        y - np.min(y),
        data_2d,
        shading="nearest",
        cmap=cmap,
    )
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return mesh


def _plot_soil_layers(ds, x, y, field_name, title_prefix, output_path, nlevels=4):
    if field_name not in ds:
        logger.info("Skipping %s plot: variable not present", field_name)
        return

    arr = np.asarray(ds[field_name].values)
    if arr.ndim != 3:
        logger.warning(
            "Skipping %s plot: expected 3D (z,y,x), got shape %s",
            field_name,
            arr.shape,
        )
        return

    n_layers = min(nlevels, arr.shape[0])
    ncols = 2
    nrows = int(math.ceil(n_layers / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows), squeeze=False)
    flat_axes = axes.ravel()

    for level in range(n_layers):
        mesh = _plot_2d_field(
            flat_axes[level],
            x,
            y,
            arr[level, :, :],
            f"{title_prefix} level {level + 1}",
        )
        fig.colorbar(mesh, ax=flat_axes[level])

    for idx in range(n_layers, len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.tight_layout()
    fig.savefig(output_path / f"{field_name}_levels.png", dpi=300)
    plt.close(fig)


def _get_cover_weighted_variable_names(ds, short_names):
    """Return sorted base variable names that can be cover-weighted.

    A variable is considered cover-weightable when for at least one short
    name there is both ``cover_<short>`` and ``<base>_<short>`` in the dataset.
    """

    cover_weighted = set()
    data_var_names = list(ds.data_vars.keys())

    for var_name in data_var_names:
        for short in short_names:
            suffix = f"_{short}"
            if not var_name.endswith(suffix):
                continue

            base_name = var_name[: -len(suffix)]
            if base_name == "cover":
                continue

            cover_name = f"cover_{short}"
            if cover_name in ds and var_name in ds:
                cover_weighted.add(base_name)

    return sorted(cover_weighted)


def _plot_all_cover_weighted_spatial_maps(ds, x, y, short_names, output_path):
    """Plot one spatial map per base variable using cover-weighted blending.

    For each base variable ``var``, computes:
        var_weighted(x, y) = sum_i(cover_i * var_i) / sum_i(cover_i)
    where ``i`` spans all cover short names that provide ``var_i``.
    """

    weighted_vars = _get_cover_weighted_variable_names(ds, short_names)
    if not weighted_vars:
        logger.info("No cover-weighted LSM variables found to plot")
        return

    weighted_output = output_path / "lsm_weighted_vars"
    os.makedirs(weighted_output, exist_ok=True)

    for base_name in weighted_vars:
        weighted_sum = None
        cover_sum = None

        for short in short_names:
            cover_name = f"cover_{short}"
            var_name = f"{base_name}_{short}"
            if cover_name not in ds or var_name not in ds:
                continue

            cover = np.asarray(ds[cover_name].values, dtype=float)
            values = np.asarray(ds[var_name].values, dtype=float)

            if cover.ndim != 2 or values.ndim != 2:
                logger.info(
                    "Skipping contribution %s for %s due to non-2D shapes cover=%s var=%s",
                    short,
                    base_name,
                    cover.shape,
                    values.shape,
                )
                continue

            if weighted_sum is None:
                weighted_sum = np.zeros_like(values, dtype=float)
                cover_sum = np.zeros_like(cover, dtype=float)

            weighted_sum += cover * values
            cover_sum += cover

        if weighted_sum is None or cover_sum is None:
            continue

        with np.errstate(divide="ignore", invalid="ignore"):
            blended = np.where(cover_sum > 0, weighted_sum / cover_sum, np.nan)

        fig, ax = plt.subplots(figsize=(8, 6))
        mesh = _plot_2d_field(ax, x, y, blended, f"Cover-weighted {base_name}")
        fig.colorbar(mesh, ax=ax)
        fig.tight_layout()
        fig.savefig(weighted_output / f"{base_name}.png", dpi=300)
        plt.close(fig)


def _plot_lcz_variables(lsm_netcdf_path, output_path):
    lsm_name = pathlib.Path(lsm_netcdf_path).name
    if not lsm_name.startswith("lsm.inp_"):
        logger.info("Could not infer LCZ dataset from %s", lsm_name)
        return

    lcz_name = lsm_name.replace("lsm.inp_", "lcz_ds_")
    lcz_path = pathlib.Path(lsm_netcdf_path).parent / lcz_name

    if not lcz_path.exists():
        logger.info("Skipping LCZ plots: file does not exist at %s", lcz_path)
        return

    ds_lcz = xr.open_dataset(lcz_path, engine="netcdf4")
    os.makedirs(output_path / "lcz", exist_ok=True)

    x = ds_lcz["x"].values if "x" in ds_lcz.coords else None
    y = ds_lcz["y"].values if "y" in ds_lcz.coords else None

    for var_name, da in ds_lcz.data_vars.items():
        values = np.asarray(da.values)

        if values.ndim == 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            if x is not None and y is not None and values.shape == (y.size, x.size):
                mesh = _plot_2d_field(ax, x, y, values, f"LCZ: {var_name}")
            else:
                mesh = ax.imshow(values, origin="lower", cmap="viridis")
                ax.set_title(f"LCZ: {var_name}")
            fig.colorbar(mesh, ax=ax)
            fig.tight_layout()
            fig.savefig(output_path / "lcz" / f"{var_name}.png", dpi=300)
            plt.close(fig)
            continue

        if values.ndim == 1:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(values)
            ax.set_title(f"LCZ: {var_name}")
            ax.set_xlabel("index")
            ax.set_ylabel(var_name)
            ax.grid(alpha=0.3)
            fig.tight_layout()
            fig.savefig(output_path / "lcz" / f"{var_name}.png", dpi=300)
            plt.close(fig)
            continue

        if values.ndim == 3:
            n_slices = values.shape[0]
            ncols = min(3, n_slices)
            nrows = int(math.ceil(n_slices / ncols))
            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(5 * ncols, 4 * nrows),
                squeeze=False,
            )
            flat_axes = axes.ravel()

            for idx in range(n_slices):
                arr = values[idx, :, :]
                ax = flat_axes[idx]
                if x is not None and y is not None and arr.shape == (y.size, x.size):
                    mesh = _plot_2d_field(ax, x, y, arr, f"{var_name}[{idx}]")
                else:
                    mesh = ax.imshow(arr, origin="lower", cmap="viridis")
                    ax.set_title(f"{var_name}[{idx}]")
                fig.colorbar(mesh, ax=ax)

            for idx in range(n_slices, len(flat_axes)):
                flat_axes[idx].axis("off")

            fig.tight_layout()
            fig.savefig(output_path / "lcz" / f"{var_name}.png", dpi=300)
            plt.close(fig)
            continue

        logger.info(
            "Skipping LCZ variable %s with unsupported shape %s",
            var_name,
            values.shape,
        )

    ds_lcz.close()


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

    ds = xr.open_dataset(lsm_netcdf_path, engine="netcdf4")
    cover_names = []
    short_names = []
    cover_vars = []

    for luname, lushort in zip(ds["luname"].values, ds["lushort"].values):
        short_name = _decode_text(lushort)
        long_name = _decode_text(luname)
        cover_val = ds[f"cover_{short_name}"]
        if np.sum(cover_val) > 0:
            cover_names.append(long_name)
            short_names.append(short_name)
            cover_vars.append(f"cover_{short_name}")

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
    fig.savefig(plot_base_path / "cover.png", dpi=300)
    plt.close()

    # Soil fields: 4 subplots (or fewer if dataset has fewer levels)
    _plot_soil_layers(
        ds,
        x,
        y,
        field_name="t_soil",
        title_prefix="Soil temperature",
        output_path=pathlib.Path(plot_base_path),
        nlevels=4,
    )
    _plot_soil_layers(
        ds,
        x,
        y,
        field_name="theta_soil",
        title_prefix="Soil moisture",
        output_path=pathlib.Path(plot_base_path),
        nlevels=4,
    )

    # Plot one cover-weighted spatial map for each LSM variable family
    # (e.g. ar, br, lai, lambda_s, z0m, etc.) from per-cover fields.
    _plot_all_cover_weighted_spatial_maps(
        ds,
        x,
        y,
        short_names=short_names,
        output_path=pathlib.Path(plot_base_path),
    )

    # Plot every variable available in the matching LCZ output dataset
    _plot_lcz_variables(
        lsm_netcdf_path=lsm_netcdf_path, output_path=pathlib.Path(plot_base_path)
    )

    ds.close()
    plt.ion()
    return


if __name__ == "__main__":
    # Example usage
    example_lsm_netcdf_path = pathlib.Path(
        "/Users/andrevanginkel/Documents/40_Input_and_Runs/42_Dales_Cases/42.01_generated_cases/002_LSM_tests/input/lsm.inp_001.nc"
    )
    example_plot_base_path = pathlib.Path(
        "/Users/andrevanginkel/Documents/40_Input_and_Runs/42_Dales_Cases/42.01_generated_cases/002_LSM_tests/profiles"
    )
    os.makedirs(example_plot_base_path, exist_ok=True)
    plot_lsm_cover(example_lsm_netcdf_path, example_plot_base_path)
