class openbc_counts_idx:
    def __init__(self, grid):
        # i1, j1, k1, i2, j2, imax, jmax, kmax, nxturb=0, nyturb=0, nzturb=0
        # boundary(west)%nx1turb = nyturb; boundary(west)%nx2turb = nzturb
        boundaries = {
            "west": {
                "nx1": jmax,
                "nx2": kmax,
                "nx1u": jmax,
                "nx2u": kmax,
                "nx1v": j1,
                "nx2v": kmax,
                "nx1w": jmax,
                "nx2w": k1,
                "nx1turb": nyturb,
                "nx2turb": nzturb,
            },
            # boundary(east)%nx1turb = nyturb; boundary(2)%nx2turb = nzturb
            "east": {
                "nx1": jmax,
                "nx2": kmax,
                "nx1u": jmax,
                "nx2u": kmax,
                "nx1v": j1,
                "nx2v": kmax,
                "nx1w": jmax,
                "nx2w": k1,
                "nx1turb": nyturb,
                "nx2turb": nzturb,
            },
            # boundary(south)%nx1turb = nxturb; boundary(3)%nx2turb = nzturb
            "south": {
                "nx1": imax,
                "nx2": kmax,
                "nx1u": i1,
                "nx2u": kmax,
                "nx1v": imax,
                "nx2v": kmax,
                "nx1w": imax,
                "nx2w": k1,
                "nx1turb": nxturb,
                "nx2turb": nzturb,
            },
            # boundary(north)%nx1turb = nxturb; boundary(4)%nx2turb = nzturb
            "north": {
                "nx1": imax,
                "nx2": kmax,
                "nx1u": i1,
                "nx2u": kmax,
                "nx1v": imax,
                "nx2v": kmax,
                "nx1w": imax,
                "nx2w": k1,
                "nx1turb": nxturb,
                "nx2turb": nzturb,
            },
            # boundary(top)%nx1turb = nxturb; boundary(5)%nx2turb = nyturb
            "top": {
                "nx1": imax,
                "nx2": jmax,
                "nx1u": i1,
                "nx2u": jmax,
                "nx1v": imax,
                "nx2v": j1,
                "nx1w": imax,
                "nx2w": jmax,
                "nx1turb": nxturb,
                "nx2turb": nyturb,
            },
        }
        # u0(2:i2,2:j1,1:kmax)
        # v0(2:i1,2:j2,1:kmax)
        # w0(2:i1,2:j1,1:k1)
        # thl0(2:i1,2:j1,1:kmax)
        # qt0(2:i1,2:j1,1:kmax)
        # e120(2:i1,2:j1,1:kmax)
        # sv0(2:i1,2:j1,1:kmax,tracer_prop(n)%trac_idx)
        openbc_initfields_idx = {
            "u0": {"x": (2, i2), "y": (2, j1), "z": (1, kmax)},
            "v0": {"x": (2, i1), "y": (2, j2), "z": (1, kmax)},
            "w0": {"x": (2, i1), "y": (2, j1), "z": (1, k1)},
            "thl0": {"x": (2, i1), "y": (2, j1), "z": (1, kmax)},
            "qt0": {"x": (2, i1), "y": (2, j1), "z": (1, kmax)},
            "e120": {"x": (2, i1), "y": (2, j1), "z": (1, kmax)},
            "sv0": {"x": (2, i1), "y": (2, j1), "z": (1, kmax)},
        }

        var_initfield_counts = {
            var: {
                dim: (
                    openbc_initfields_idx[var][dim][1]
                    - openbc_initfields_idx[var][dim][0]
                    + 1
                )
                for dim in ["x", "y", "z"]
            }
            for var in openbc_initfields_idx
        }

        var_boundary_counts = {
            boundary: {
                "u0": (boundaries[boundary]["nx1u"], boundaries[boundary["nx2u"]]),
                "v0": (boundaries[boundary]["nx1v"], boundaries[boundary["nx2v"]]),
                "w0": (boundaries[boundary]["nx1w"], boundaries[boundary["nx2w"]]),
                "thl0": (boundaries[boundary]["nx1"], boundaries[boundary["nx2"]]),
                "qt0": (boundaries[boundary]["nx1"], boundaries[boundary["nx2"]]),
                "e120": (boundaries[boundary]["nx1"], boundaries[boundary["nx2"]]),
                "sv0": (boundaries[boundary]["nx1"], boundaries[boundary["nx2"]]),
                "u2": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "v2": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "w2": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "uv": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "uw": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "vw": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "thl2": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "qt2": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "wthl": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
                "wqt": (
                    boundaries[boundary]["nx1turb"],
                    boundaries[boundary["nx2turb"]],
                ),
            }
            for boundary in boundaries
        }
