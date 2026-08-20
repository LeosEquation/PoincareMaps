@testset "PoincareMap - harmonic oscillator" begin

    # Harmonic oscillator
    #
    # q' = p
    # p' = -q
    #
    # Hamiltonian:
    #
    # H(q,p) = (q² + p²)/2

    function f!(dx, x, params, t)
        dx[1] = x[2]
        dx[2] = -x[1]
    end

    function H(x, params, t)
        return 0.5 * (x[1]^2 + x[2]^2)
    end

    # No coordinate transformation is required.
    function bc!(x, params, t)
        nothing
    end

    # Poincaré section:
    #
    # q = 0
    #
    # Only crossings with p > 0 are accepted.
    function g(dx, x, params, t)
        return (x[2] > 0, x[1])
    end

    # PoincaréMap searches the second coordinate (p)
    # for the value satisfying H(q,p) = E.
    function nrfind(x, params, t)
        idx = 2

        p_min = 0.0
        p_max = 1.0

        Δp = (p_max - p_min) / 50

        p_min += Δp
        p_max -= Δp

        return idx, Δp, p_max, p_min
    end

    # Two initial seeds with different q values.
    seed = [
        [-0.5, 0.0],
        [0.5, 0.0],
    ]

    E = 0.5
    params = nothing

    t0 = 0.0
    tmax = 100.0

    abstol = 1e-16
    order = 20

    crossings = PoincareMap(
        f!,
        bc!,
        g,
        H,
        nrfind,
        seed,
        E,
        params,
        tmax,
        abstol,
        order;
        maxsteps=1000,
        maxevents=500,
        eventorder=0,
        newtoniter=20,
        nrabstol=1e-14,
        Etol=1e-14,
    )

    # At least one crossing must have been detected.
    @test size(crossings, 1) > 0

    # The Poincaré section is q = 0.
    q_error = maximum(abs.(crossings[:, 1]))

    @test q_error < 1e-12

    # At E = 1/2, every crossing with p > 0 satisfies
    #
    # H(0,p) = p²/2 = 1/2
    #
    # and therefore p = 1.
    p_error = maximum(abs.(crossings[:, 2] .- 1.0))

    @test p_error < 1e-12

end