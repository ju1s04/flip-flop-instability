# Bondi–Hoyle–Lyttleton Accretion in Castro

**Problem directory:** `Castro-main/Exec/science/bhl_accretion/`
**Goal:** implement 2-D Bondi–Hoyle–Lyttleton (BHL) / Hoyle–Lyttleton (HL) wind accretion
onto a point mass in Castro, targeting a replication of the **flip-flop instability**
study of **Blondin & Pope (2009)**, in particular their **Model A**.

---

## 1. Physical problem

A uniform, supersonic gas wind (density `ρ∞`, speed `V∞`, sound speed `c∞`,
Mach number `ℳ = V∞/c∞`) flows past a gravitating point mass `M`. Gravity
focuses the flow into a bow shock and a downstream accretion wake/column. The
characteristic scale is the **Hoyle–Lyttleton accretion radius**

```
Ra = 2 G M / V∞²                                   (Hoyle & Lyttleton 1939)
```

the impact parameter inside which incoming gas is captured. The analytic
accretion rates are

```
Ṁ_HL = π Ra² ρ∞ V∞           (3-D / cylindrical, Hoyle & Lyttleton 1939)
Ṁ_2D = 2  Ra  ρ∞ V∞          (2-D planar, Blondin & Pope 2009, eq. 8)
```

In 2-D planar geometry the accretion wake is subject to the **flip-flop
instability** (Fryxell & Taam 1988): the wake swings from side to side with an
exponentially growing amplitude, periodically forming transient counter-rotating
accretion disks. Blondin & Pope (2009) characterise this as a true *overstability*
and fit the accreted specific angular momentum `j(t)` with an exponentially
growing sinusoid,

```
j(t) = j0 · exp(ω_r t) · cos(ω_i t)                (their eq. 9)
```

extracting a growth rate `ω_r` and oscillation frequency `ω_i`. Their **Model A**
(the canonical case) gives `ω_r = 0.070`, `ω_i = 0.68`, period `2π/ω_i = 9.2 ≈
2π√2 = 8.9` (the Keplerian period at `Ra`), and a steady pre-instability accretion
rate `Ṁ ≈ 0.9 Ṁ_2D`.

---

## 2. Model parameters (this implementation)

We work in **dimensionless units** with `ρ∞ = 1`, `V∞ = 1`, and `Ra = 1`
(so `GM = 0.5`). Time is measured in units of `Ra / V∞ = 1`.

| Quantity | Symbol | Value (code units) | Notes |
|---|---|---|---|
| Ambient density | `ρ∞` | **1.0** | |
| Wind speed | `V∞` | **1.0** | +x direction |
| Mach number | `ℳ` | **4.0** | Blondin & Pope Model A |
| Sound speed | `c∞ = V∞/ℳ` | **0.25** | |
| Adiabatic index | `γ` | **4/3 = 1.3333** | ideal gas |
| Ambient pressure | `p∞ = ρ∞c∞²/γ` | **0.046875** | |
| Gravitational parameter | `GM` | **0.5** | |
| Accretion radius | `Ra = 2GM/V∞²` | **1.0** | length unit |
| Castro `point_mass` | `M = GM/G` | **7.49143×10⁶** | `G=6.67428×10⁻⁸` (CGS), so `G·M = 0.5` |
| 2-D accretion rate | `Ṁ_2D = 2Raρ∞V∞` | **2.0** | validation target `0.9×2 = 1.8` |
| 3-D accretion rate | `Ṁ_HL = πRa²ρ∞V∞` | **π ≈ 3.1416** | |

### Accretor / sink size

- **Accretor (sink) radius:** `Rs = Rs_frac × Ra`, with `Rs_frac = 0.25`
  → **`Rs = 0.25 Ra`** in the current default.
- **Blondin & Pope Model A uses `Rs = 0.037 Ra`** (a *small* accretor). Reaching
  that value is the outstanding task (see §8): it requires both a working
  absorbing sink and `dx ≪ 0.037 Ra` (≈ 4–5 AMR levels).
- **Important:** all figures produced so far were run with the **sink disabled**
  (`castro.add_ext_src = 0`). In those runs there is *no resolved accretor
  surface*,  only the gravitating point at the domain centre. The black marker in
  the zoom figure is the point-mass location, not a resolved sink.

### Domain and grid

- **Domain:** `[-6, 18] × [-8, 8]` in units of `Ra` (24 Ra × 16 Ra). Accretor at
  the centre, `(x, y) = (6, 0)`. Upstream extent 12 Ra, downstream 12 Ra,
  transverse ±8 Ra.
- **Base grid (default):** `768 × 512` → `dx = 0.03125 Ra`. (Fast test runs used a
  `192 × 128` base → `dx = 0.125 Ra`.)
- **AMR:** `max_level = 2` default; production zoom used 3 levels. Finest cell
  size reached so far: **`dx = 0.0156 Ra`** (192^2 base + 3 levels).

---

## 3. Castro implementation

The problem is under work so it can be self-contained on CASTRO's `Exec/science/` directory. Note that **no core Castro
source was modified**. Files:

| File | Role |
|---|---|
| `GNUmakefile` | `DIM=2`, `USE_MPI=TRUE`, `USE_GRAV=TRUE`, `USE_REACT=FALSE`, `EOS_DIR=gamma_law`, `NETWORK_DIR=general_null` (single passive species) |
| `_prob_params` | runtime parameters (see §2) |
| `problem_initialize.H` | sets `problem::center` to the domain centre; derives `Ra`, `Rs`, `c∞`, `p∞`, `ρe∞`, `T∞` |
| `problem_initialize_state_data.H` | **initial condition**: uniform Mach-4 wind everywhere |
| `problem_bc_fill.H` | **custom inflow BC**: imposes the uniform wind on the XLO face |
| `problem_source.H` | the **absorbing sink** (floor-relaxation drain inside `Rs`) — *currently disabled* |
| `problem_tagging.H` | **AMR tagging**: nested, distance-based refinement around the accretor |
| `inputs` | full runtime configuration |
| `plot_bhl.py`, `plot_amr.py`, `plot_zoom.py` | yt/matplotlib post-processing |

### 3.1 Numerical methods

| Aspect | Choice | Castro control |
|---|---|---|
| Integrator | **CTU** (corner transport upwind) + Strang | `time_integration_method = 0` (default) |
| Reconstruction | **PPM** (Colella & Woodward 1984) | `ppm_type = 1` |
| Riemann solver | **CGF** (Colella–Glaz–Ferguson, two-shock) | `riemann_solver = 0` |
| EOS | ideal gas | `gamma_law`, `eos_gamma = 1.3333`, `eos_assume_neutral = 1` |
| Energy | dual-energy formalism (tracks `ρe` separately) | Castro default |
| Robustness | timestep retry; flux limiting near low density | `use_retry = 1`, `limit_fluxes_on_small_dens = 1` |
| CFL | 0.5 | `cfl = 0.5` |

### 3.2 Boundary conditions

`castro.lo_bc = 1 2`, `castro.hi_bc = 2 2` (codes: `1 = Inflow`, `2 = Outflow`).

- **XLO = Inflow (upstream):** Dirichlet, the uniform Mach-4 wind is imposed via
  the custom `problem_bc_fill.H` (modelled on `hydro_tests/double_mach_reflection`).
- **XHI = Outflow (downstream):** zero-gradient; the supersonic wake exits here.
- **YLO / YHI = Outflow:** transverse far field.


### 3.3 Gravity

**Native point-mass gravity** (no self-gravity, no Poisson solve):

```
castro.do_grav        = 1
gravity.gravity_type  = ConstantGrav
gravity.const_grav    = 0.0
castro.use_point_mass = 1
castro.point_mass     = 7.49143e6
```

Castro applies `g = −G · point_mass·r̂ / r²` toward `problem::center`, coupled into
the CTU predictor/corrector with a gravity-limited timestep. This configuration
mirrors Castro's own point-mass accretion problem `hydro_tests/rotating_torus`.
The gravity is **un-softened** (raw `1/r²`); the accretor sits on a cell corner so
`r ≥ dx/√2 > 0`, but `g` grows as `dx` shrinks on finer AMR levels.

### 3.4 The absorbing sink (status: written, disabled)

`problem_source.H` implements a floor-relaxation absorbing sink: inside `Rs` it
drains the *excess* of each conserved quantity above a positive floor
(`ρ_sink = rho_sink_frac × ρ∞ = 10⁻³`), with a smooth taper and a drain rate capped
at `0.5/dt`. It is enabled by `castro.add_ext_src = 1`. **It is currently disabled
(`= 0`)** because it destabilises near the un-softened gravitational cusp
(see §6). Mirroring Blondin & Pope, the intended fix is a *hard reset* of the
absorbed region rather than a fractional drain (§8).

### 3.5 AMR strategy

Castro has no box-based refinement in the inputs file, so `problem_tagging.H`
tags cells by **distance from the accretor**, with a radius that **halves per
level** (nested refinement):

```
r_lim(level) = refine_radius / 2^level
```

With `refine_radius = 3`: level 1 covers `r < 3 Ra`, level 2 `r < 1.5 Ra`, level 3
`r < 0.75 Ra`, concentrating the finest cells tightly on the accretor. This both
resolves the accretion flow and keeps the finest level affordable.

### 3.6 Provenance of key design choices (where each choice came from)

For reproducibility, this records *why* the two least-obvious choices — the
gravity treatment and the boundary conditions — were made, and their sources.

**Gravity = native point-mass (`use_point_mass`).** Three sources, in order:

1. *Physics — why a point mass at all:* **Blondin & Pope (2009), §2**, which
   states the model uses *"a point source of gravity at the coordinate origin."*
   BHL accretion is onto an external compact object with negligible gas
   self-gravity, so a point mass (not Poisson self-gravity) is the correct model.
2. *Implementation — the exact Castro config:* copied from Castro's own
   point-mass accretion test, **`Exec/hydro_tests/rotating_torus/inputs_3d`**,
   which uses `gravity.gravity_type = ConstantGrav`, `const_grav = 0.0`,
   `use_point_mass = 1`. This choice was reached *after* an initial attempt to
   apply softened gravity through `problem_source.H` (a generic `ext_src`) crashed
   at the bow shock (§6, item 5); searching for Castro's own working accretion
   setup led to `rotating_torus`.
3. *Units — the value `point_mass = 7.49143×10⁶`:* read directly from the kernel
   **`Source/gravity/Gravity.cpp`** (~line 3028),
   `radial_force = -C::Gconst * castro::point_mass / rsq`, with
   `Gconst = 6.67428×10⁻⁸` (CGS, hard-coded). Setting `point_mass = GM/Gconst`
   with `GM = 0.5` gives the intended dimensionless field.

> Note: Castro *does* have a full Poisson self-gravity solver (`PoissonGrav`,
> multigrid/Hypre) and a `MonopoleGrav` approximation. They are deliberately
> **not** used here — self-gravity of the accreting gas is negligible in BHL and
> is not part of the HL model, and it would add the expensive `solve_for_phi`
> Poisson kernel for no physical benefit.

**Boundary conditions (XLO inflow, others outflow).** Also three sources:

1. *Physics:* wind-past-an-object flow ⇒ Dirichlet inflow upstream, zero-gradient
   outflow downstream and transverse. Same logic as Blondin & Pope (inflow at the
   upstream/outer boundary, outflow downstream), on a Cartesian grid.
2. *Castro conventions:* the BC integer codes
   (`0=Interior, 1=Inflow, 2=Outflow, 3=Symmetry, 4=SlipWall`) and the
   custom-inflow mechanism were taken from existing problems —
   `hydro_tests/RT/inputs_2d` (codes) and `hydro_tests/double_mach_reflection`
   (the spatially-varying `problem_bc_fill.H` template).
3. *A constraint learned from a crash:* attempting two Inflow faces meeting at a
   corner aborted with *"external boundaries meeting at a corner not supported,"*
   which is why only the upstream face is Inflow. A pure-wind isolation run
   (gravity + sink off) then confirmed the BCs were correct and that the earlier
   crash was the gravity source, not the boundaries.

---

## 4. Simulations performed & results

All runs use `mpirun -np 4` on the local machine (macOS, Apple Silicon, open-mpi +
homebrew GCC 15). Native point-mass gravity, **sink disabled**.

| Run | Grid | Result | Figure |
|---|---|---|---|
| Uniform | `384²`-class, `max_level = 0` | **stable** to `t = 1.63 Ra/V`, 17 plotfiles | `bhl_density.png` |
| AMR ×2 | `192×128` base + 2 levels | **stable** to `t = 1.8 Ra/V` | `bhl_amr.png` |
| AMR ×3 | `192×128` base + 3 levels, `dx=0.0156 Ra` | **stable** to `t = 1.23 Ra/V` | `bhl_zoom.png` |

**Physics captured (qualitative):** the uniform Mach-4 wind develops a bow shock,
a gravitationally-focused dense accretion structure, and a downstream wake. By
`t ≈ 1.2–1.8` the wake is **asymmetric** and shows a rotating spiral/vortex — the
**onset of the flip-flop instability** (transient rotating accretion structure).
Peak density reaches ~30–77× ambient at the accretor (the pile-up is expected
because the sink is off — nothing removes the accreted gas).

**Figures**
- `bhl_density.png` — uniform-grid time series (`t = 0.19, 0.60, 1.20, 1.63`).
- `bhl_amr.png` — full domain + zoom with AMR grid boxes (levels 0/1/2).
- `bhl_zoom.png` — 2×2 Ra close-up at the finest level with velocity streamlines,
  showing the bow shock and the flip-flop vortex.

---

## 5. Local build recipe (macOS, Apple Silicon)

`mpicxx` wraps Apple clang by default; point open-mpi at homebrew GCC so
`COMP=gnu` gets a matched toolchain, and **build twice** the first time
(GNU Make 3.81 races the generated `runtime_params.cpp`):

```bash
export OMPI_CC=gcc-15 OMPI_CXX=g++-15 OMPI_FC=gfortran
make COMP=gnu USE_MPI=TRUE -j8      # pass 1: generates sources (may error on the race)
make COMP=gnu USE_MPI=TRUE -j8      # pass 2: links the executable
mpirun -np 4 ./Castro2d.gnu.MPI.ex inputs
```

---


## 6. Comparison with Blondin & Pope (2009)

| Paper element (their §2–3) | This work |
|---|---|
| 2-D HL accretion, point mass, ideal gas |  same model |
| Mach 4, γ = 4/3, ε = 0 (Model A parameters) |  used |
| Source terms: **gravity only** |  gravity only (sink off = no absorption) |
| **PPM** reconstruction (Colella & Woodward) |  PPM (`ppm_type=1`) |
| Dual-energy / internal-energy handling (their §2.1) |  Castro's built-in dual energy |
| Grid: 2-D **polar (r,φ)** | **Cartesian + AMR** (Castro has no planar-polar; see below) |
| **Absorbing inner boundary** at `Rs` | not working (sink disabled) |
| **Ballistic** upstream BC (their §2.2) | uniform-wind inflow |
| **Steady-state** initial condition (their §2.3) | uniform-flow start |
| **Angular-momentum-conserving** scheme (their §2.5) | Castro conserves linear momentum |
| Qualitative flip-flop onset (their Figs 1, 4, 8) | captured (asymmetric wake + spiral) |
| **Steady `Ṁ ≈ 0.9 Ṁ_2D`** | not yet measured |
| **Growth rate `ω_r`, Table 1** | not yet measured |

**On polar coordinates.** AMReX supports only `cartesian`, `RZ` (axisymmetric
r–z) and `SPHERICAL` (axisymmetric r–θ). Blondin & Pope use **planar polar
(r,φ)**, which Castro does not have, and which would require adding a coordinate
type to AMReX plus polar geometric source terms throughout Castro's hydro; a codebase development still under work. 


---

## 7. Current status, limitations, and roadmap to Model A

**Working now:** build (local MPI), initialisation, boundary conditions, native
point-mass gravity, AMR (nested tagging, ≥3 levels), plotting, all stable, and
the flip-flop *onset* is visible.

**Not yet done:** any *quantitative* Blondin & Pope result. The two gating pieces:

1. **A robust absorbing accretor.** Replace the fractional-drain sink with a
   **hard reset** inside `R_core` (pin `ρ, v, p` to fixed absorbed values every
   step.
2. **The `j(t)` diagnostic**, accreted specific angular momentum, normalised to
   `Rs Vc = √(GM Rs)` — to fit `ω_r`.

---

## References

- Beckmann, R. S., et al. 2018, *Bondi or not Bondi: resolution and accretion/drag
  for supermassive black holes*, MNRAS 478, 995.
- Benensohn, J. S., Lamb, D. Q., & Taam, R. E. 1997, ApJ 478, 723.
- Blondin, J. M., & Pope, T. C. 2009, *Revisiting the "Flip-Flop" Instability of
  Hoyle-Lyttleton Accretion*, ApJ 700, 95 (arXiv:0905.2769). **[primary reference]**
- Blondin, J. M., & Raymer, E. 2012, *Hoyle-Lyttleton Accretion in Three
  Dimensions*, ApJ 752, 30.
- Bondi, H., & Hoyle, F. 1944, MNRAS 104, 273.
- Colella, P., & Woodward, P. R. 1984, *The Piecewise Parabolic Method (PPM)*,
  J. Comput. Phys. 54, 174.
- Edgar, R. 2004, *A Review of Bondi–Hoyle–Lyttleton Accretion*,
  New Astron. Rev. 48, 843.
- Foglizzo, T., Galletti, P., & Ruffert, M. 2005, A&A 435, 397.
- Fryxell, B. A., & Taam, R. E. 1988, ApJ 335, 862 (FT — flip-flop discovery).
- Hoyle, F., & Lyttleton, R. A. 1939, Proc. Camb. Phil. Soc. 35, 405 & 592.
- Hunt, R. 1971, MNRAS 154, 141 (axisymmetric HL, `Ṁ = 0.88 Ṁ_HL`).
- Krumholz, M. R., McKee, C. F., & Klein, R. I. 2004, ApJ 611, 399 (sink particles).
- Matsuda, T., Inoue, M., & Sawada, K. 1987, MNRAS 226, 785 (flip-flop flow).
- Ruffert, M., & Arnett, D. 1994, ApJ 427, 351; Ruffert, M. 1999, A&A 346, 861
  (3-D nested-Cartesian HL accretion).
- Shima, E., et al. 1985, MNRAS 217, 367; 1998, A&A 337, 311.
- Almgren, A. S., et al. 2010, *CASTRO: A New Compressible Astrophysical Solver*,
  ApJ 715, 1221.
- Zhang, W., et al. 2019, *AMReX*, J. Open Source Softw. 4, 1370.

---

*Report generated for the `bhl_accretion` Castro problem. Figures: `bhl_density.png`,
`bhl_amr.png`, `bhl_zoom.png`. Post-processing: `plot_bhl.py`, `plot_amr.py`,
`plot_zoom.py`.*
