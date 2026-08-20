# src/poincaremap.jl
function PoincareMap(
    f!,
    bc!,
    g,
    H,
    nrfind,
    seed::Array{Vector{U}},
    E::U,
    params,
    tmax::U,
    abstol::U,
    order::T;
    parse_eqs::Bool=true,
    maxsteps::Int=500,
    maxevents::Int=500,
    eventorder::Int=0,
    newtoniter::Int=10,
    nrabstol::U=eps(U),
    Etol::U=eps(E),
) where {U<:Real,T<:Int}

    dim = length(seed[1])
    N = length(seed)

    # progress = Threads.Atomic{Int}(0)
    p = 0
    step = ceil(Int, N / 9)

    t0 = zero(U)

    ξ = seed[1]
    idx, Δξ, ξ_max, ξ_min = nrfind(ξ, params, t0)
    ξ[idx] = ξ_min
    mindim = 1
    (mindim <= idx <= dim) || error("Índice inválido: $idx")
    (ξ_min < ξ_max) || error("Valores máximo y mínimo no válidos.")
    (Δξ <= (ξ_max - ξ_min)) || error("Paso no válido")

    caches = [init_cache_ps(t0, ξ, maxevents, order, f!, params; parse_eqs) for _ in 1:N]
    ξs = Array{U}(undef, dim, N)
    ξauxs = [Taylor1(0.0, 1) for _ in 1:dim, _ in 1:N]
    Hauxs = [Taylor1(0.0, 1) for _ in 1:N]

    for i in 1:N
        for j in 1:dim
            ξs[j, i] = seed[i][j]
        end
        ξauxs[idx, i][1] = 1.0
    end

    crosses = [Array{U}(undef, dim, maxevents + 1) for _ in 1:N]

    nevents = zeros(Int, N)

    print("_Progress_\n")

    for i in 1:N

        @views ξ = ξs[:, i]
        @views ξaux = ξauxs[:, i]
        Haux = Hauxs[i]
        cache = caches[i]
        cross = crosses[i]

        idx, Δξ, ξ_max, ξ_min = nrfind(ξ, params, t0)
        (mindim <= idx <= dim) || error("Índice inválido: $idx")
        (ξ_min < ξ_max) || error("Valores máximo y mínimo no válidos.")

        ξ_old = ξ_min
        ξ_new = ξ_old + Δξ
        ξ[idx] = ξ_old

        if ξ_old < ξ_max

            E_old = H(ξ, params, t0)

            while ξ_new <= ξ_max

                ξ[idx] = ξ_new

                E_new = H(ξ, params, t0)

                if (E_old - E) * (E_new - E) < 0.0

                    for j in 1:dim
                        ξaux[j][0] = ξ[j]
                    end

                    niter = 1

                    Haux .= H(ξaux, params, t0)

                    while niter <= newtoniter &&
                        abs(Haux[0] - E) > Etol

                        δξ = (Haux[0] - E) / Haux[1]
                        ξaux[idx][0] -= δξ

                        val = ξaux[idx][0]

                        if ξ_min < val < ξ_max
                            Haux .= H(ξaux, params, t0)
                            niter += 1
                        end

                    end

                    ξ[idx] = ξaux[idx][0]

                    if abs(H(ξ, params, t0) - E) <= Etol

                        nev = taylorinteg_ps!(
                            f!, bc!, g, ξ, t0, tmax, abstol,
                            cache, params; maxsteps, maxevents, eventorder, newtoniter, nrabstol
                        ) - 1

                        if nev > 0
                            nevents[i] = nev
                            for comp1 in 1:dim
                                for comp2 in 1:nev
                                    cross[comp1, comp2] = cache.xv[comp1, comp2]
                                end
                            end
                        end

                    end
                end

                ξ_old = ξ_new
                E_old = E_new
                ξ_new += Δξ
            end

        end

        p += 1 # Threads.atomic_add!(progress, 1)

        if p % step == 0
            print("█")
        end
    end

    effective_crosses = Array{U,2}(undef, dim, sum(nevents))

    comp1 = 0
    for comp2 in 1:N
        for comp3 in 1:nevents[comp2]
            comp1 += 1
            for k in 1:dim
                effective_crosses[k, comp1] = crosses[comp2][k, comp3]
            end
        end
    end

    print("█\n")

    print("‾‾‾‾‾‾‾‾‾‾\n")

    return transpose(effective_crosses)

end