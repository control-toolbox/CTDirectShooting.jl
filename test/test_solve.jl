function test_solve()

    @testset "basic solving" begin

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

        sol = solve(ocp,30)

        @test sol.objective ≈ 2.0 atol=1e-3

    end

    @testset "basic solving 2" begin

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
        
        sol = solve(ocp,30)

        println(sol)

        #@test sol.objective ≈ 2.0 atol=1e-3
    
    end
    
end