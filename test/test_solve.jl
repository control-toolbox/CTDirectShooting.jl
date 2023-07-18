function test_solve()

    @testset "basic solving" begin
        t0 = 0
        x0 = [-1, 0]
        xf = [0, 0]
        γ = 1

        g(t0, x0, tf, xf) = tf
        f⁰(t, x, u) = 0
        f(x, u) = [x[2];u]
        dmin = [-Inf, -Inf, -γ]
        dmax = [Inf, Inf, γ]
        cmin = []
        c(t, x, u) = nothing
        cmax = []
        bmin = [0, -1, 0, 0, 0, 0]
        b(t0, x0, tf, xf) = [t0, x0, tf, xf]
        bmax = [0, -1, 0, Inf, 0, 0]

        prob = SimpleProblem(g, f⁰, f, dmin, dmax, cmin, c, cmax, bmin, b, bmax, 2, 1)

        solve(prob, 9, 3)
    end
    
end