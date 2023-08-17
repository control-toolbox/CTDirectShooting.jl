module CTDirectShooting

using LinearAlgebra
using CTFlows
using ADNLPModels
using NLPModelsIpopt
using OrdinaryDiffEq
using NLPModels

include("utils.jl")
include("nlp_solver.jl")
include("init.jl")
include("problem.jl")
include("solution.jl")
include("check.jl")
include("solve.jl")

export SimpleProblem, solve

end