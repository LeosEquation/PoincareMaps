# src/integration/taylorinteg.jl

import TaylorIntegration: update_cache!, taylorstep!

function taylorinteg_ps!(
    f!,
    bc!,
    g,
    q0::AbstractVector{U},
    t0::T,
    tmax::T,
    abstol::T,
    cache::VectorCachePS,
    params;
    maxsteps::Int=500,
    maxevents::Int=500,
    eventorder::Int=0,
    newtoniter::Int=10,
    nrabstol::T=eps(T),
    reltol::T=eps(T),
) where {T<:Real,U<:Number}

    (; xv, xaux, t, x, dx, rv, parse_eqs) = cache

    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    sign_tstep = copysign(1, tmax - t0)

    g_tupl = g(dx, x, params, t)
    g_tupl_old = deepcopy(g_tupl)
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = evaluate(x_dx)
    g_dg_val = vcat(evaluate(g_tupl[2]), evaluate(g_tupl_old[2]))

    nsteps = 1
    nevents = 1

    while sign_tstep * t0 < sign_tstep * tmax

        δt_old = δt

        δt = taylorstep!(
            Val(parse_eqs),
            f!,
            t,
            x,
            dx,
            xaux,
            abstol,
            params,
            rv,
            reltol,
        )

        if iszero(δt)
            break
        end

        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))

        evaluate!(x, δt, x0)

        g_tupl = g(dx, x, params, t)

        nevents = findroot_ps!(
            bc!,
            params,
            t0,
            x,
            dx,
            g_tupl_old,
            g_tupl,
            eventorder,
            xv,
            δt_old,
            x_dx,
            x_dx_val,
            g_dg,
            g_dg_val,
            nrabstol,
            newtoniter,
            nevents,
        )

        g_tupl_old = deepcopy(g_tupl)
        t0 += δt

        bc!(x0, params, t0)
        update_cache!(cache, t0, x0)

        nsteps += 1

        if nsteps > maxsteps || nevents > maxevents
            break
        end
    end

    return nevents
end

function taylorinteg_ps!(
    f!,
    bc!,
    g,
    lims,
    q0::AbstractVector{U},
    t0::T,
    tmax::T,
    abstol::T,
    cache::VectorCachePS,
    params;
    maxsteps::Int=500,
    maxevents::Int=500,
    eventorder::Int=0,
    newtoniter::Int=10,
    nrabstol::T=eps(T),
    reltol::T=zero(T),
) where {T<:Real,U<:Number}

    (; xv, xaux, t, x, dx, rv, parse_eqs) = cache

    x0 = deepcopy(q0)
    update_cache!(cache, t0, x0)
    sign_tstep = copysign(1, tmax - t0)

    g_tupl = g(dx, x, params, t)
    g_tupl_old = deepcopy(g_tupl)
    δt = zero(x[1])
    δt_old = zero(x[1])

    x_dx = vcat(x, dx)
    g_dg = vcat(g_tupl[2], g_tupl_old[2])
    x_dx_val = evaluate(x_dx)
    g_dg_val = vcat(evaluate(g_tupl[2]), evaluate(g_tupl_old[2]))

    nsteps = 1
    nevents = 1

    while sign_tstep * t0 < sign_tstep * tmax

        δt_old = δt

        δt = taylorstep!(
            Val(parse_eqs),
            f!,
            t,
            x,
            dx,
            xaux,
            abstol,
            params,
            rv,
            reltol,
        )

        if iszero(δt)
            break
        end

        δt = sign_tstep * min(δt, sign_tstep * (tmax - t0))

        evaluate!(x, δt, x0)

        g_tupl = g(dx, x, params, t)

        nevents = findroot_ps!(
            bc!,
            params,
            t0,
            x,
            dx,
            g_tupl_old,
            g_tupl,
            eventorder,
            xv,
            δt_old,
            x_dx,
            x_dx_val,
            g_dg,
            g_dg_val,
            nrabstol,
            newtoniter,
            nevents,
        )

        g_tupl_old = deepcopy(g_tupl)
        t0 += δt

        bc!(x0, params, t0)

        if lims(x0, params, t0)
            break
        end

        update_cache!(cache, t0, x0)

        nsteps += 1

        if nsteps > maxsteps || nevents > maxevents
            break
        end
    end

    return nevents
end