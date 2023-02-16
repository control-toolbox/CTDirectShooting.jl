using CTDirectShooting
using Test
using CTBase
using CTProblemLibrary
using ControlToolboxTools
using Plots

# CTDirectShooting
const CTOptimizationInit = CTDirectShooting.CTOptimizationInit
const convert_init = CTDirectShooting.convert_init
const __init = CTDirectShooting.__init
const __grid = CTDirectShooting.__grid

# CTBase
const __grid_size_direct_shooting = CTBase.__grid_size_direct_shooting
const __init_interpolation = CTBase.__init_interpolation

#
@testset verbose = true showtiming = true "Direct shooting" begin
    for name in (
        "CTOptimization", # unconstrained direct simple shooting
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end
