@testset "PoincareMap - nonlinear Hamiltonian" begin

    # ------------------------------------------------------------
    # Nonlinear Hamiltonian:
    #
    # H = 1/2 (p1^2 + p2^2 + q1^2 + q2^2)
    #     + α q1^2 q2^2
    #
    # State:
    # x = [q1, q2, p1, p2]
    # ------------------------------------------------------------

    α = 0.1

    function f!(dx, x, params, t)

        q1, q2, p1, p2 = x

        dx[1] = p1
        dx[2] = p2

        dx[3] = -q1 - 2α*q1*q2^2
        dx[4] = -q2 - 2α*q2*q1^2

    end


    # Hamiltonian
    function H(x, params, t)

        q1, q2, p1, p2 = x

        return (
            0.5*(p1^2 + p2^2 + q1^2 + q2^2)
            +
            α*q1^2*q2^2
        )

    end


    # No modification of the state.
    function bc!(x, params, t)
        nothing
    end


    # ------------------------------------------------------------
    # Poincaré section:
    #
    # q1 = 0
    #
    # Only crossings with p1 > 0 are accepted.
    # ------------------------------------------------------------

    function g(dx, x, params, t)

        p1 = dx[3][0]

        return (p1 > 0, x[1])

    end

    # ------------------------------------------------------------
    # Find the coordinate used to impose H = E.
    #
    # We vary p2.
    # ------------------------------------------------------------

    function nrfind(x, params, t)

        idx = 4

        Δp2 = 0.02

        p2_min = -1.5
        p2_max = 1.5

        return idx, Δp2, p2_max, p2_min

    end


    # ------------------------------------------------------------
    # Initial seeds.
    #
    # The values of p2 will subsequently be adjusted by
    # PoincareMap so that H = E.
    # ------------------------------------------------------------

    seed = [
        [0.2, 0.3, 0.1, 0.0],
        [0.4, 0.2, 0.1, 0.0],
    ]


    E = 1.0
    params = nothing

    tmax = 100.0
    abstol = 1e-14
    order = 20

    maxsteps = 5000
    maxevents = 100


    # ------------------------------------------------------------
    # Compute Poincaré map.
    # ------------------------------------------------------------

    crossings = PoincareMaps.PoincareMap(
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
        maxsteps=maxsteps,
        maxevents=maxevents,
        eventorder=0,
        newtoniter=20,
        nrabstol=1e-14,
        Etol=1e-12,
    )


    # ------------------------------------------------------------
    # Basic sanity check.
    # ------------------------------------------------------------

    @test size(crossings, 2) == 4
    @test size(crossings, 1) > 0


    # ------------------------------------------------------------
    # Every returned point must belong to the energy surface.
    # ------------------------------------------------------------

    energies = [
        H(view(crossings, i, :), params, 0.0)
        for i in axes(crossings, 1)
    ]

    @test maximum(abs.(energies .- E)) < 1e-10


    # ------------------------------------------------------------
    # Every returned point must belong to the Poincaré section:
    #
    # q1 ≈ 0
    # ------------------------------------------------------------

    q1_error = maximum(abs.(crossings[:, 1]))

    @test q1_error < 1e-10


    # ------------------------------------------------------------
    # Every returned crossing must satisfy the condition imposed
    # by g:
    #
    # p1 > 0
    # ------------------------------------------------------------

    @test all(crossings[:, 3] .> 0)

end