import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.patches import Patch
import logging

logger = logging.getLogger(__name__)
logger.debug("Entered module: %s", __name__)

ds = xr.open_dataset(
    "/Users/andre/Documents/Documenten/PhD/40_Input_and_Runs/42_Dales_Cases/Case_generator/02_32x32_grass_with_lake/input/lsm.inp_001.nc"
)
# Example: categories and colors
categories = ["cover_fbd", "cover_grs", "cover_urb", "cover_aqu"]
colors = {
    "cover_fbd": "tab:green",  # e.g. forest
    "cover_grs": "tab:olive",  # grass
    "cover_urb": "tab:gray",  # urban
    "cover_aqu": "tab:blue",  # water
}

# Extract coordinate grids
x = ds["x"].values
y = ds["y"].values
nx, ny = len(x), len(y)

dx = x[1] - x[0]
dy = y[1] - y[0]

fig, ax = plt.subplots(figsize=(5, 5))

# Loop through each grid cell
for j, yval in enumerate(y):
    for i, xval in enumerate(x):
        bottom = 0
        for cat in categories:
            frac = float(ds[cat].isel(x=i, y=j))
            rect = plt.Rectangle(
                (xval - 0.5 * dx, yval - 0.5 * dy + bottom),  # lower-left corner
                dx,  # cell width
                frac * dy,  # fraction height (fills vertically)
                facecolor=colors[cat],
                edgecolor="none",
            )
            ax.add_patch(rect)
            bottom += frac * dy

# Formatting
ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
ax.set_ylim(y.min() - 0.5, y.max() + 0.5)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")

legend_patches = [
    Patch(facecolor=color, label=cat.replace("_frac", ""))
    for cat, color in colors.items()
]
ax.legend(handles=legend_patches, loc="upper right")
# ax.legend(frameon=False, loc='upper right')

plt.tight_layout()
plt.savefig("test.png")
plt.close()
