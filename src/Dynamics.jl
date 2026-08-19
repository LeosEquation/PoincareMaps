using TaylorIntegration, TaylorSeries, LinearAlgebra,
    PeriodicSchurDecompositions, Printf

# Inicializando objetos
include("./objects/equilibrium_bifurcations.jl")
include("./objects/equilibrium_branch.jl")
include("./objects/limitpoint_branch.jl")
include("./objects/branchpoint_branch.jl")
include("./objects/orbit_bifurcations.jl")
include("./objects/orbit_branch.jl")
include("./objects/flippoint_branch.jl")
include("./inicialization/taylornini.jl")

# include("./ploting/ploteqbranch.jl")

# Inicializando operaciones extra
# include("./operations/adjugate.jl")
# include("./operations/biproduct.jl")
# include("./operations/modangle.jl")

# Equilibrio
include("./equilibria/eqfinding.jl")
include("./equilibria/eqsystem.jl")
include("./equilibria/eqcontinuation.jl")
include("./equilibria/eqbranchswitching.jl")

# Limit Points
include("./limitpoint/lpfinding.jl")
include("./limitpoint/lpsystem.jl")
include("./limitpoint/lpcontinuation.jl")
# include("./limitpoint/hlpcontinuation.jl")

# # Hopf
# include("./hopf/hpffinding.jl")
# include("./hopf/hpfsystem.jl")
# include("./hopf/hpfcontinuation.jl")

# Branch Points
include("./branchpoint/bpfinding.jl")
include("./branchpoint/bpsystem.jl")
include("./branchpoint/bpcontinuation.jl")
# include("./branchpoint/hbpcontinuation.jl")

# Orbits 
# include("./orbits/orbitfinding.jl")
include("./orbits/auxiliary.jl")
include("./orbits/orbitsystem.jl")
include("./orbits/orbitcontinuation.jl")
include("./orbits/orbitcontinuation_extend.jl")
include("./orbits/dporbitcontinuation.jl")
# include("./orbits/multishooting.jl")
# include("./orbits/horbitcontinuation.jl")
# include("./orbits/horbitcontinuation_initorbit.jl")
# include("./orbits/hbranchswitching.jl")

# Branch Cycle Point
include("./branchpointcycle/bpcfinding.jl")
include("./branchpointcycle/bpcsystem.jl")
include("./branchpointcycle/bpccontinuation.jl")

# Flip Point
include("./flippoint/fpsystem.jl")
include("./flippoint/fpfinding.jl")
include("./flippoint/fpcontinuation.jl")
include("./flippoint/fpcontinuation_extend.jl")
# include("./perioddoubling/pdsystem.jl")
# include("./perioddoubling/pdcontinuation.jl")

# Poincare Section 
include("./poincaresection/poincaresection.jl")
# include("./poincaresection/energylims.jl")
include("./poincaresection/energyfinding.jl")