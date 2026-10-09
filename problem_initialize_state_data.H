#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    amrex::ignore_unused(i, j, k, geomdata);

    // Uniform supersonic wind moving in the +x direction.
    // (This is the simplest initial condition; see the README for the
    //  steady-state / ballistic-inflow refinements used by Blondin & Pope.)

    state(i,j,k,URHO)  = problem::rho_inf;
    state(i,j,k,UMX)   = problem::rho_inf * problem::v_inf;
    state(i,j,k,UMY)   = 0.0_rt;
    state(i,j,k,UMZ)   = 0.0_rt;
    state(i,j,k,UEINT) = problem::rhoe_inf;
    state(i,j,k,UEDEN) = problem::rhoe_inf
                       + 0.5_rt * problem::rho_inf * problem::v_inf * problem::v_inf;
    state(i,j,k,UTEMP) = problem::T_inf;

    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = problem::rho_inf * problem::xn[n];
    }
}

#endif
