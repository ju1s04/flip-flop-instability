#!/usr/bin/env python3
"""Tight zoom on the accretor at the finest AMR level, with streamlines."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt
yt.set_log_level(50)

FN = "plt00300"           # most-evolved frame
X0, X1 = 5.0, 7.0         # tighter zoom box (Ra), accretor at (6,0)
Y0, Y1 = -1.0, 1.0

ds = yt.load(FN)
lev = ds.index.max_level
dxf = ds.domain_width.d[0] / (ds.domain_dimensions[0] * ds.refine_by**lev)
nx = int(round((X1 - X0) / dxf))
ny = int(round((Y1 - Y0) / dxf))

# covering grid over ONLY the zoom box, at the finest level -> full L2 detail
le = [X0, Y0, float(ds.domain_left_edge.d[2])]
cg = ds.covering_grid(level=lev, left_edge=le, dims=[nx, ny, 1])
rho = np.array(cg["density"][:, :, 0]).T
vx  = np.array((cg["xmom"][:, :, 0] / cg["density"][:, :, 0])).T
vy  = np.array((cg["ymom"][:, :, 0] / cg["density"][:, :, 0])).T

extent = [X0, X1, Y0, Y1]
fig, ax = plt.subplots(figsize=(9.5, 8), constrained_layout=True)
im = ax.imshow(np.log10(rho), origin="lower", extent=extent, aspect="equal",
               cmap="magma", vmin=-0.3, vmax=1.5)

xs = np.linspace(X0, X1, nx)
ys = np.linspace(Y0, Y1, ny)
ax.streamplot(xs, ys, vx, vy, color="cyan", density=1.6,
              linewidth=0.6, arrowsize=0.8)
ax.plot(6.0, 0.0, marker="o", ms=11, mfc="black", mec="white", mew=1.0)
ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
ax.set_xlabel("x / Ra"); ax.set_ylabel("y / Ra")
ax.set_title(f"BHL accretion, finest level (dx={dxf:.4f} Ra), "
             f"t = {float(ds.current_time):.2f} Ra/V")
plt.colorbar(im, ax=ax, shrink=0.85, label=r"$\log_{10}\rho$")
fig.savefig("bhl_zoom.png", dpi=140)
print(f"wrote bhl_zoom.png  ({nx}x{ny} at level {lev}, dx={dxf:.4f})")
print(f"zoom rho range: [{rho.min():.3f}, {rho.max():.3f}]")
