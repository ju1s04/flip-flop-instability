#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // place the accretor at the center of the domain
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }

    // single passive species (pure ideal gas)
    for (int n = 0; n < NumSpec; ++n) {
        problem::xn[n] = 0.0_rt;
    }
    problem::xn[0] = 1.0_rt;

    // derived ambient (upstream) state
    problem::c_inf = problem::v_inf / problem::mach;
    problem::p_inf = problem::rho_inf * problem::c_inf * problem::c_inf / problem::gamma_gas;

    // Hoyle-Lyttleton accretion radius:  Ra = 2 G M / V_inf^2
    problem::Ra  = 2.0_rt * problem::GM / (problem::v_inf * problem::v_inf);
    problem::Rs  = problem::Rs_frac * problem::Ra;
    problem::eps = problem::soft_frac * problem::Rs;

    // specific internal energy and temperature of the ambient gas
    eos_t eos_state;
    eos_state.rho = problem::rho_inf;
    eos_state.p   = problem::p_inf;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = problem::xn[n];
    }
    eos(eos_input_rp, eos_state);

    problem::rhoe_inf = problem::rho_inf * eos_state.e;
    problem::T_inf    = eos_state.T;

    amrex::Print() << "BHL accretion: Ra = " << problem::Ra
                   << ", Rs = " << problem::Rs
                   << ", Mdot_2D = 2 Ra rho_inf v_inf = "
                   << 2.0_rt * problem::Ra * problem::rho_inf * problem::v_inf
                   << std::endl;
}

#endif
