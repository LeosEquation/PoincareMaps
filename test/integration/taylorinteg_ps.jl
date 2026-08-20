# test/integration/taylorinteg_ps.jl

# test/integration/taylorinteg_ps.jl

@testset "taylorinteg_ps! - harmonic oscillator" begin

    # Harmonic oscillator:
    #
    # q' = p
    # p' = -q
    #
    # Exact solution for q(0)=1, p(0)=0:
    #
    # q(t) = cos(t)
    # p(t) = -sin(t)

    function f!(dx, x, params, t)
        dx[1] = x[2]
        dx[2] = -x[1]
    end

    # Boundary condition: no modification of the state.
    function bc!(x, params, t)
        nothing
    end

    # Poincaré section:
    #
    # g(x) = q
    #
    # The second component must be a Taylor polynomial because
    # findroot_ps! operates on its Taylor coefficients.
    function g(dx, x, params, t)
        return (true, x[1])
    end

    q0 = [1.0, 0.0]
    params = nothing

    t0 = 0.0
    tmax = 10π

    abstol = 1e-16
    order = 20

    maxevents = 20
    maxsteps = 1000

    cache = PoincareMaps.init_cache_ps(
        t0,
        q0,
        maxevents,
        order,
        f!,
        params;
        parse_eqs=true,
    )

    nevents = PoincareMaps.taylorinteg_ps!(
        f!,
        bc!,
        g,
        q0,
        t0,
        tmax,
        abstol,
        cache,
        params;
        maxsteps=maxsteps,
        maxevents=maxevents,
        eventorder=0,
    )

    # The implementation starts nevents at 1, so the number
    # of actual detected crossings is nevents - 1.
    ncrossings = nevents - 1

    @test ncrossings == 10

    crossings = cache.xv[:, 1:ncrossings]

    q_error = maximum(abs.(crossings[1, :]))

    @test q_error < 1e-12

    p_expected = [
        isodd(i) ? -1.0 : 1.0
        for i in 1:ncrossings
    ]

    p_error = maximum(
        abs.(crossings[2, :] .- p_expected)
    )

    @test p_error < 1e-12

end