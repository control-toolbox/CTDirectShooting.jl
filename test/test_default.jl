
    @testset "Default value of the grid size for the direct shooting method" begin
        @test CTBase.__grid_size_direct_shooting() isa Integer
    end

    @testset "Default value of the penalty term for the direct shooting method" begin
        @test CTBase.__penalty_term_direct_shooting() isa Real
    end

    @testset "Default value of the maximum number of iterations for the direct shooting method" begin
        @test CTBase.__max_iter_direct_shooting() isa Integer
    end

    @testset "Default value of the absolute tolerance for the direct shooting method" begin
        @test CTBase.__abs_tol_direct_shooting() isa Real
    end

    @testset "Default value of the optimality tolerance for the direct shooting method" begin
        @test CTBase.__opt_tol_direct_shooting() isa Real
    end

    @testset "Default value of the stagnation tolerance for the direct shooting method" begin
        @test CTBase.__stagnation_tol_direct_shooting() isa Real
    end
