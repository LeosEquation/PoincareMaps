# src/integration/cache.jl

import TaylorIntegration: AbstractVectorCache, init_expansions, _determine_parsing!, update_cache!

struct VectorCachePS{XV,XAUX,T,X,DX,RV,PARSE_EQS} <: AbstractVectorCache
    xv::XV
    xaux::XAUX
    t::T
    x::X
    dx::DX
    rv::RV
    parse_eqs::PARSE_EQS
end

function init_cache_ps(
    t0::T,
    q0::AbstractVector{U},
    maxevents::Int,
    order::Int,
    f!,
    params=nothing;
    parse_eqs::Bool=true,
) where {U,T}
    # Initialize the vector of Taylor1 expansions
    t, x, dx = init_expansions(t0, q0, order)
    # Determine if specialized jetcoeffs! method exists
    parse_eqs, rv = _determine_parsing!(parse_eqs, f!, t, x, dx, params)
    # Initialize cache
    dof = length(q0)
    return VectorCachePS(
        Array{U}(undef, dof, maxevents + 1),
        Array{Taylor1{U}}(undef, dof),
        t,
        x,
        dx,
        rv,
        parse_eqs,
    )
end

function update_cache!(cache::AbstractVectorCache, t0::T, x0::AbstractVector{U}) where {T,U}
    (; t, x) = cache
    @inbounds for i in eachindex(x0)
        x[i][0] = x0[i]
    end
    @inbounds t[0] = t0
    return nothing
end
