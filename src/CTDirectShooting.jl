module CTDirectShooting

#
using CTBase
using CTFlows
using CTOptimization
using LinearAlgebra # for the norm for instance
using MLStyle
using Printf # to print iterations results for instance

# Other declarations
const nlp_constraints = CTBase.nlp_constraints
const __grid_size_direct_shooting = CTBase.__grid_size_direct_shooting
const __display = CTBase.__display
const __penalty_constraint = CTBase.__penalty_constraint
const __callbacks = CTBase.__callbacks
const __init_interpolation = CTBase.__init_interpolation
const __iterations = CTBase.__iterations 
const __absoluteTolerance = CTBase.__absoluteTolerance
const __optimalityTolerance = CTBase.__optimalityTolerance
const __stagnationTolerance = CTBase.__stagnationTolerance
const expand = CTBase.expand
const vec2vec = CTBase.vec2vec

#
include("utils.jl")
include("init.jl")
include("problem.jl")
include("solution.jl")
include("solve.jl")

# export functions only for user
export solve

end