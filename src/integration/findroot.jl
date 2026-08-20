# src/integration/findroot.jl

import TaylorIntegration: surfacecrossing, nrconvergencecriterion

function findroot_ps!(
    bc!,
    params,
    t0,
    x,
    dx,
    g_tupl_old,
    g_tupl,
    eventorder,
    xvS,
    δt_old,
    x_dx,
    x_dx_val,
    g_dg,
    g_dg_val,
    nrabstol,
    newtoniter,
    nevents,
)

    if surfacecrossing(g_tupl_old, g_tupl, eventorder)
        #auxiliary variables
        g_val = g_tupl[2]
        g_val_old = g_tupl_old[2]
        nriter = 1
        dof = length(x)

        #first guess: linear interpolation
        slope = (g_val[eventorder] - g_val_old[eventorder]) / δt_old
        dt_li = -(g_val[eventorder] / slope)

        x_dx[1:dof] = x
        x_dx[(dof+1):2dof] = dx
        g_dg[1] = derivative(g_val, eventorder)
        g_dg[2] = derivative(g_dg[1])

        #Newton-Raphson iterations
        dt_nr = dt_li
        evaluate!(g_dg, dt_nr, view(g_dg_val, :))

        while nrconvergencecriterion(g_dg_val[1], nrabstol, nriter, newtoniter)
            dt_nr = dt_nr - g_dg_val[1] / g_dg_val[2]
            evaluate!(g_dg, dt_nr, view(g_dg_val, :))
            nriter += 1
        end

        if nriter <= newtoniter

            evaluate!(x_dx, dt_nr, view(x_dx_val, :))

            bc!(x_dx_val, params, t0 + dt_nr)

            xvS[:, nevents] .= deepcopy.(view(x_dx_val, 1:dof))

            nevents += 1

        end
    end

    return nevents
end