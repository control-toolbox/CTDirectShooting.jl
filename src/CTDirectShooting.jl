"""
[`CTDirectShooting`](@ref) module.

Lists all the imported modules and packages:

$(IMPORTS)

List of all the exported names:

$(EXPORTS)

"""
module CTDirectShooting

using LinearAlgebra
using CTFlows
using ADNLPModels
using NLPModelsIpopt
using OrdinaryDiffEq
using NLPModels
using DocStringExtensions
using CTBase
using StaticArrays

include("utils.jl")
include("nlp_solver.jl")
include("init.jl")
include("problem.jl")
include("solution.jl")
include("check.jl")
include("solve.jl")

export SimpleProblem, solve

end