#!/usr/bin/env python3
"""Plot the AMR BHL run: full-domain density + a zoom with AMR grid boxes."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import yt
yt.set_log_level(50)

ds = yt.load("plt00300")
lev = ds.index.max_level
cg = ds.covering_grid(level=lev, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions * ds.refine_by**lev)
rho = np.array(cg["density"][:, :, 0]).T
le = ds.domain_left_edge.d
re = ds.domain_right_edge.d
extent = [le[0], re[0], le[1], re[1]]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)

im1 = ax1.imshow(np.log10(rho), origin="lower", extent=extent, aspect="equal",
                 cmap="magma", vmin=-1, vmax=1.5)
ax1.plot(6, 0, "c+", ms=12, mew=2)
ax1.set_title(f"log density  (t = {float(ds.current_time):.2f} Ra/V)")
ax1.set_xlabel("x / Ra"); ax1.set_ylabel("y / Ra")
plt.colorbar(im1, ax=ax1, shrink=0.8, label=r"$\log_{10}\rho$")

im2 = ax2.imshow(np.log10(rho), origin="lower", extent=extent, aspect="equal",
                 cmap="magma", vmin=-1, vmax=1.5)
colors = ["white", "cyan", "yellow", "red"]
for g in ds.index.grids:
    L = int(g.Level)
    gle = g.LeftEdge.d; gre = g.RightEdge.d
    ax2.add_patch(Rectangle((gle[0], gle[1]), gre[0]-gle[0], gre[1]-gle[1],
                            fill=False, edgecolor=colors[L % 4], lw=0.7))
ax2.plot(6, 0, "c+", ms=12, mew=2)
ax2.set_xlim(3.5, 9.5); ax2.set_ylim(-3, 3)
ax2.set_title("zoom + AMR grids (white=L0, cyan=L1, yellow=L2)")
ax2.set_xlabel("x / Ra"); ax2.set_ylabel("y / Ra")
plt.colorbar(im2, ax=ax2, shrink=0.8, label=r"$\log_{10}\rho$")

fig.suptitle("Bondi-Hoyle-Lyttleton accretion with AMR (native point-mass gravity)",
             fontsize=13)
fig.savefig("bhl_amr.png", dpi=130)
print("wrote bhl_amr.png")
print("AMR levels present:", sorted(set(int(g.Level) for g in ds.index.grids)))
print("cells per level:", {L: sum(g.ActiveDimensions.prod() for g in ds.index.grids if int(g.Level)==L)
                           for L in sorted(set(int(g.Level) for g in ds.index.grids))})
