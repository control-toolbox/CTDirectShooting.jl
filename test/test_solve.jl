function test_solve()

    @testset "basic solving" begin
        t0 = 0.0
        x0 = [-1.0, 0.0]
        xf = [0.0, 0.0]
        γ = 1.0

        g(t0, x0, tf, xf) = tf
        f⁰(t, x, u) = 0.0
        f(t, x, u) = [x[2];u]
        # function f(t,x::AbstractVector{T},u) where {T <: Real}
        #     ẋ = Vector{T}(undef, size(x,1))
        #     println("                                                  ",typeof(ẋ))
        #     ẋ[1] = x[2]
        #     ẋ[2] = u
        #     return ẋ
        # end
        dmin = [-Inf, -Inf, -γ]
        dmax = [Inf, Inf, γ]
        cmin = [0.0]
        c(t, x, u) = [0.0]
        cmax = [0.0]
        bmin = [0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
        b(t0, x0, tf, xf) = [t0; x0; tf; xf]
        bmax = [0.0, -1.0, 0.0, Inf, 0.0, 0.0]

        prob = SimpleProblem(g, f⁰, f, dmin, dmax, cmin, c, cmax, bmin, b, bmax, 2, 1)

        sol = solve(prob, 30, 3)

        println(sol.solution)
    end
    
end