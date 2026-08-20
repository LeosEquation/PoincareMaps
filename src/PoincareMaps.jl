# src/PoincareMaps.jl

module PoincareMaps

using TaylorIntegration
using TaylorSeries
using LinearAlgebra

# Integration
include("integration/cache.jl")
include("integration/findroot.jl")
include("integration/taylorinteg.jl")

# Poincaré map
include("poincaremap.jl")

export PoincareMap

end