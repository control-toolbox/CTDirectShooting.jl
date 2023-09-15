function test_solve()

    @testset "basic solving : double integrator final time" begin

        t0=0
        x0=[-1, 0]
        xf=[0, 0]
        γ = 1
        A = [ 0 1
            0 0 ]
        B = [ 0
        1 ]

        @def ocp begin
            tf ∈ R, variable
            t ∈ [ t0, tf ], time
            x ∈ R², state
            u ∈ R, control
            x(t0) == x0,    (initial_con) 
            x(tf) == xf,    (final_con)
            -γ ≤ u(t) ≤ γ,  (u_con)
            ẋ(t) == A * x(t) + B * u(t)
            tf → min
        end

        sol = solve(ocp,grid_size_fine=30)

        @test sol.objective ≈ 2.0 atol=1e-3

    end

    @testset "basic solving : simple exponential final time" begin

        t0 = 0
        x0 = 0

        @def ocp begin
            tf ∈ R, variable
            t ∈ [t0, tf], time
            x ∈ R, state
            u ∈ R, control
            x(t0) == x0, initial_con
            (x(tf) - tf) - 10 == 0, boundary_constraint
            ẋ(t) == u(t)
            ∫(0.5 * u(t) ^ 2) → min
        end
        
        sol = solve(ocp,grid_size_fine=30,grid_size_coarse=10,fixed_time_step=true)

        @test sol.objective ≈ 20.0 atol=1e-3
    
    end

    @testset "basic solving : double integrator energy" begin

        t0=0
        tf=1
        x0=[-1, 0]
        xf=[0, 0]
        A = [ 0 1
            0 0 ]
        B = [ 0
        1 ]

        @def ocp begin
            t ∈ [t0, tf], time
            x ∈ R², state
            u ∈ R, control
            x(t0) == [-1, 0], initial_con
            x(tf) == [0, 0], final_con
            ẋ(t) == A * x(t) + B * u(t)
            ∫(0.5 * u(t) ^ 2) → min
        end

        sol = solve(ocp,grid_size_fine=30)

        @test sol.objective ≈ 6.0 atol=1e-2
    
    end
    
end