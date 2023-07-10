using CTDirectShooting
using Test

# CTDirectShooting


#
@testset verbose = true showtiming = true "Direct shooting" begin
    for name in (
         # unconstrained direct simple shooting
        )
        @testset "$name" begin
            include("test_$name.jl")
        end
    end
end
