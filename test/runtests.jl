using CTDirectShooting
using Test

# CTDirectShooting


#
@testset verbose = true showtiming = true "Direct shooting" begin
    for name ∈ (
        :init,
        :utils,
        :solve,
        :problem,
         # unconstrained direct simple shooting
        )
        @testset "$(name)" begin
            test_name = Symbol(:test_, name)
            include("$(test_name).jl")
            @eval $test_name()
        end
    end
end
