#ifndef problem_tagging_H
#define problem_tagging_H

#include <prob_parameters.H>

// Refine a fixed region around the accretor so the accretion flow (and the
// sink) are always resolved at the finest level.  Castro has no box-based
// refinement in the inputs file, so we tag by distance from problem::center.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    amrex::ignore_unused(state);
    amrex::ignore_unused(level);

    const Real* problo = geomdata.ProbLo();
    const Real* dx     = geomdata.CellSize();

    Real xx = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
    Real yy = 0.0_rt;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
    Real r = std::sqrt(xx * xx + yy * yy);

    // Level-dependent (nested) refinement: the tag radius halves with each
    // level, so the finest level hugs the accretor.  Level L tags cells for
    // refinement to L+1, so:  L1 covers refine_radius, L2 covers R/2,
    // L3 covers R/4, ...
    Real rlim = problem::refine_radius / std::pow(2.0_rt, static_cast<Real>(level));

    if (r < rlim) {
        tag(i,j,k) = TagBox::SET;
    }
}

#endif
