#ifndef problem_source_H
#define problem_source_H

#include <prob_parameters.H>

using namespace amrex;

// Absorbing accretion sink: smoothly drain gas inside Rs.
//
// Gravity is NOT applied here -- it is handled by Castro's native point-mass
// gravity (castro.use_point_mass), which couples correctly to the CTU hydro
// and supplies a gravity-limited timestep.  This source only removes mass,
// momentum and energy inside the sink radius, capped so it can never remove
// more than half a cell's content in one step.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{
    amrex::ignore_unused(time);

    const Real* problo = geomdata.ProbLo();
    const Real* dx     = geomdata.CellSize();

    Real xx = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
    Real yy = 0.0_rt;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
    Real r = std::sqrt(xx * xx + yy * yy);

    // default: no source
    src(i,j,k,URHO)  = 0.0_rt;
    src(i,j,k,UMX)   = 0.0_rt;
    src(i,j,k,UMY)   = 0.0_rt;
    src(i,j,k,UMZ)   = 0.0_rt;
    src(i,j,k,UEDEN) = 0.0_rt;
    src(i,j,k,UEINT) = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        src(i,j,k,UFS+n) = 0.0_rt;
    }

    if (r < problem::Rs) {
        Real r_core = 0.5_rt * problem::Rs;
        Real f_sink;
        if (r <= r_core) {
            f_sink = 1.0_rt;
        } else {
            Real t = (problem::Rs - r) / (problem::Rs - r_core);
            f_sink = t * t;
        }

        // cap the drain so we never remove more than ~half a cell per step
        Real inv_tau = amrex::min(1.0_rt / problem::tau_sink, 0.5_rt / dt);
        Real rate = f_sink * inv_tau;

        // Relax toward an absorbed state: a low-but-positive floor density,
        // zero velocity, and the corresponding internal energy.  Only the
        // EXCESS above the floor is removed, so density/energy can never reach
        // exactly zero (which the un-softened point-mass gravity would trigger).
        Real rho_sink  = problem::rho_sink_frac * problem::rho_inf;
        Real e_spec    = problem::rhoe_inf / problem::rho_inf;   // ambient specific internal energy
        Real rhoe_sink = rho_sink * e_spec;

        src(i,j,k,URHO)  = -amrex::max(state(i,j,k,URHO)  - rho_sink,  0.0_rt) * rate;
        src(i,j,k,UMX)   = -state(i,j,k,UMX) * rate;
        src(i,j,k,UMY)   = -state(i,j,k,UMY) * rate;
        src(i,j,k,UMZ)   = -state(i,j,k,UMZ) * rate;
        src(i,j,k,UEDEN) = -amrex::max(state(i,j,k,UEDEN) - rhoe_sink, 0.0_rt) * rate;
        src(i,j,k,UEINT) = -amrex::max(state(i,j,k,UEINT) - rhoe_sink, 0.0_rt) * rate;
        for (int n = 0; n < NumSpec; ++n) {
            src(i,j,k,UFS+n) = -amrex::max(state(i,j,k,UFS+n) - rho_sink, 0.0_rt) * rate;
        }
    }
}

#endif
