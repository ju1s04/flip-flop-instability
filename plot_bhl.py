#!/usr/bin/env python3
"""Plotting the uniform-grid BHL accretion flow (density)"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

yt.set_log_level(50) 

frames = ["plt00050", "plt00150", "plt00300", "plt00400"]

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

for ax, fn in zip(axes.flat, frames):
    ds = yt.load(fn)
    lev = ds.index.max_level
    dims = (ds.domain_dimensions * ds.refine_by**lev)
    cg = ds.covering_grid(level=lev, left_edge=ds.domain_left_edge, dims=dims)
    rho = np.array(cg["density"][:, :, 0]).T   # (ny, nx)

    le = ds.domain_left_edge.d
    re = ds.domain_right_edge.d
    extent = [le[0], re[0], le[1], re[1]]

    im = ax.imshow(np.log10(rho), origin="lower", extent=extent,
                   aspect="equal", cmap="magma", vmin=-1, vmax=1.5)
    ax.plot(6.0, 0.0, "c+", ms=10, mew=2)         # accretor
    ax.set_title(f"{fn}   t = {float(ds.current_time):.2f} Ra/V", fontsize=10)
    ax.set_xlabel("x / Ra"); ax.set_ylabel("y / Ra")
    plt.colorbar(im, ax=ax, label=r"$\log_{10}\rho$", shrink=0.8)

fig.suptitle("Bondi-Hoyle-Lyttleton accretion (uniform grid, Mach 4) log density",
             fontsize=13)
fig.savefig("bhl_density.png", dpi=130)
