#ifndef problem_bc_fill_H
#define problem_bc_fill_H

// Custom boundary fill: impose the uniform supersonic wind on every inflow
// (ext_dir) face.  The setup uses inflow on XLO (upstream) and on YLO/YHI
// (transverse far field, which is just undisturbed wind), and outflow on XHI
// (downstream, where the accretion wake leaves the grid).  Imposing the wind
// on all three inflow faces -- including the shared corners -- avoids the
// inflow/outflow corner ambiguity that can leave a ghost cell with garbage.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
    amrex::ignore_unused(time);

    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();

    bool fill = false;

    if (bcs(URHO).lo(0) == amrex::BCType::ext_dir && i < domlo[0]) {
        fill = true;                       // XLO (upstream)
    }
#if AMREX_SPACEDIM >= 2
    if (bcs(URHO).lo(1) == amrex::BCType::ext_dir && j < domlo[1]) {
        fill = true;                       // YLO (transverse)
    }
    if (bcs(URHO).hi(1) == amrex::BCType::ext_dir && j > domhi[1]) {
        fill = true;                       // YHI (transverse)
    }
#endif

    if (fill) {
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
}

#endif
