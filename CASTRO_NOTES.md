# `bhl_accretion`

Two-dimensional **Bondi–Hoyle–Lyttleton (BHL) accretion**: a uniform supersonic
wind flows past a gravitating point mass equipped with an absorbing sink. This
setup reproduces the classic Hoyle–Lyttleton accretion flow and, for a
sufficiently small accretor, the 2D **"flip-flop" instability** of the accretion
wake.

## Replication target

This problem is built to reproduce:

* **Blondin & Pope (2009), ApJ, 700, 95**,
  *"Revisiting the Flip-Flop Instability of Hoyle–Lyttleton Accretion"*
  (arXiv:0905.2769)

with the analytic accretion rate and review from

* **Edgar (2004), New Astron. Rev., 48, 843** — BHL theory and Ṁ formulae,
* **Hoyle & Lyttleton (1939); Bondi & Hoyle (1944)** — original theory,
* **Fryxell & Taam (1988)** — discovery of the 2D flip-flop instability.

The flip-flop is an **inherently two-dimensional** overstability: it is absent in
3D simulations (Ruffert 1999; Blondin & Raymer 2012), so this is a genuinely-2D
problem, not a 1D flow on a 2D grid.

## Model and units

The flow is parameterized by the accretor mass `M`, wind speed `V∞`, ambient
density `ρ∞`, ambient sound speed `c∞` (i.e. Mach number `ℳ = V∞/c∞`) and the
accretor radius `Rs`. The Hoyle–Lyttleton accretion radius is

```
Ra = 2 G M / V∞²
```

The `inputs` file uses dimensionless units `ρ∞ = 1`, `V∞ = 1`, `Ra = 1`
(so `GM = 0.5`), with time measured in units of `Ra/V∞`.

## Implementation (what this problem adds)

Everything is done with existing Castro hooks — no `Source/` changes:

* `problem_initialize_state_data.H` — uniform wind initial condition.
* `problem_source.H` — **softened point-mass gravity** `g = −GM r/(r²+ε²)^{3/2}`
  plus a **smooth absorbing sink** that drains gas inside `Rs`
  (enabled by `castro.add_ext_src = 1`; `USE_GRAV = FALSE`).
* `problem_bc_fill.H` — the upstream wind **inflow** boundary on `XLO`
  (modeled on `hydro_tests/double_mach_reflection`).

## Build & run

```bash
make -j
mpirun -np 4 ./Castro2d.gnu.MPI.ex inputs
```

### Local macOS build (Apple Silicon, open-mpi + homebrew GCC)

`mpicxx` wraps Apple clang by default; point it at homebrew GCC so `COMP=gnu`
gets a matched toolchain, and build twice the first time (GNU Make 3.81 races
the generated `runtime_params.cpp`):

```bash
export OMPI_CC=gcc-15 OMPI_CXX=g++-15 OMPI_FC=gfortran
make COMP=gnu USE_MPI=TRUE -j8      # first pass generates sources
make COMP=gnu USE_MPI=TRUE -j8      # second pass links the exe
```

## Current status (verified locally)

* **Build + init + BCs + I/O**: working (MPI, 2D).
* **Native point-mass gravity + AMR**: **verified stable.** The uniform-grid and
  2-level-AMR runs both complete cleanly to t ~ 1.8 Ra/V and reproduce the
  gravitationally-focused accretion flow (bow shock + wake, becoming asymmetric
  — the flip-flop precursor). See `plot_bhl.py` / `plot_amr.py`.
* **Absorbing sink** (`castro.add_ext_src=1`): **work in progress.** The sink
  drives near-accretor cells to *exactly* zero density (which the un-softened
  point-mass gravity aggravates), tripping `enforce_min_density`. The default
  `inputs` therefore ships with the sink **off**. Options to finish it:
  (a) measure the accretion rate as the mass flux through a control surface at
  Rs instead of an absorbing sink (needs no sink at all); or
  (b) a softened-gravity + floor-relaxation sink that keeps density strictly
  positive throughout the strong-gravity region; or
  (c) Blondin & Pope's inner-boundary + steady-state-relaxation approach.

## Validation

**Level 1 — steady accretion rate.** Run the default (uniform grid, `Rs=0.25 Ra`).
The mass accretion rate should settle to

```
Ṁ ≈ 0.9 Ṁ_2D ,   Ṁ_2D = 2 Ra ρ∞ V∞
```

in agreement with Hunt (1971) and Blondin & Pope. Ṁ is the volume integral of
the sink drain rate (see "Diagnostics" below); a resolution study
(512²→768²→1024²) should show it converge.

**Level 2 — flip-flop instability (Blondin & Pope Table 1).** Switch on the AMR
block and set `problem.Rs_frac = 0.037` (model A). Measure the accreted specific
angular momentum `j(t)` (or the wake oscillation) and fit
`j(t) = j₀ e^{ω_r t} cos(ω_i t)`. Targets:

| Model | Rs/Ra | ℳ | γ | ω_r (published) |
|:-----:|:-----:|:-:|:-----:|:---------------:|
| A | 0.037  | 4  | 4/3  | 0.070 |
| B | 0.037  | 10 | 4/3  | 0.088 |
| D | 0.0125 | 4  | 4/3  | 0.104 |
| E | 0.0037 | 4  | 4/3  | 0.170 |
| G | 0.037  | 4  | 1.40 | 0.068 |
| H | 0.037  | 4  | 1.50 | 0.040 |
| I | 0.037  | 4  | 1.55 | 0.031 |
| F | 0.037  | 4  | 5/3  | 0.000 (stable) |

Also: oscillation period `2π/ω_i ≈ 9.2`, comparable to the Keplerian period at
`Ra` (`2π√2 = 8.9`); the flow is stable for `γ ≥ 1.6` at `Rs = 0.037 Ra`.

## Known differences from Blondin & Pope (honest caveats)

* **Grid geometry.** Blondin & Pope use a 2D **polar** grid with an inner radial
  absorbing boundary; this problem uses a **Cartesian** grid with a source-term
  sink. Reproducing the same period/growth on a Cartesian grid is a useful
  **grid-independence** cross-check (the flip-flop's grid dependence has been
  debated in the literature).
* **Initial condition.** They relax a half-domain **steady state** first, so the
  instability grows cleanly from small amplitude. This problem starts from
  uniform flow; the early transient (tail-shock formation) masks the very early
  growth, so for a clean `ω_r` fit you may want to relax a steady state first.
* **Upstream boundary.** They use a **ballistic** inflow (their eqs. 4–7); this
  problem uses a simpler uniform-wind inflow (accurate to ~0.5% upstream).
* **Angular momentum.** VH-1 was modified to conserve angular momentum for the
  disk phase; Castro conserves linear momentum. This affects the late disk
  dynamics, not the early growth-rate measurement.
* **Resolution.** The timestep near the sink scales as `Rs^{3/2}`, so small
  accretors are expensive — start with model A (`Rs = 0.037`) and use AMR with a
  fixed refinement region on the accretor.

## Diagnostics (TODO / next step)

Ṁ is currently obtained by post-processing (integrate `ρ · f_sink / τ` over the
sink region from plotfiles, or track the total grid mass printed by
`castro.sum_interval`). A `Problem_Derive` field exposing the local sink rate —
so Ṁ can be summed directly — is a natural next addition.
