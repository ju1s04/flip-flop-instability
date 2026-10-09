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

* **Edgar (2004), New Astron. Rev., 48, 843** BHL theory and Ṁ formulae,
* **Hoyle & Lyttleton (1939); Bondi & Hoyle (1944)** original theory,
* **Fryxell & Taam (1988)** discovery of the 2D flip-flop instability.

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

