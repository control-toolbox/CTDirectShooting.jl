function test_utils()

    @test CTDirectShooting.rg(1,1) == 1

    @testset "Traduction OCP to SimpleProblem" begin
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

        prob  = CTDirectShooting.OCP_To_SimpleProblem(ocp)
        println(prob)

    end

    
    
end