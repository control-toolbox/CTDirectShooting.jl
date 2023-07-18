using CTDirectShooting
using CTFlows
using ADNLPModels
using NLPModelsIpopt 
using OrdinaryDiffEq

function solve_example()
    # the problem (double integrator time control constraint)
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

    # solving 
    ff(x,u,t) = f(x,u)
    fv(v) = [f(v[1:2],v[3]);0]
    tf = t0+1
    N = 9
    I = [i for i ∈ 0:N]
    T = range(t0, tf, N+1)
    β = fill(γ/2,N)
    function u(t)
        index = 1
        while (t >= T[index] && index < N)
            index += 1
        end
        return(β[index])
    end
    U = [u(t) for t in T[1:end-1]]
    J = I[1:3:size(I,1)]
    M = size(J,1)
    XJ = [(((lastindex(J)-j)/size(J,1))*x0 + (j-1)/size(J,1)*xf) for j in eachindex(J)] 
    F = tf

    # variables are (tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1])  := unk

    obj_fun(unk) = unk[N+1] #
    dim_unk = size(I,1) + 2*size(J,1) + size(I,1)-1 # N+1 + M + N
    unk0 = fill(0.5,dim_unk)
    l_var = -Inf*ones(dim_unk)
    u_var = Inf*ones(dim_unk)
    for i ∈ 1:dim_unk # variable bounds
        if i ≤ N+1
            l_var[i] = 0
            u_var[i] = Inf
            #
            #unk0[i] = tf*i/(N+1)
        end
    end

    lb = zeros(6 + 2*M + 3(N+1))
    ub = zeros(6 + 2*M + 3(N+1))
    lb[1] = t0
    lb[2:3] = x0
    lb[4] = 0
    ub[5:6] = xf
    ub[1] = t0
    ub[2:3] = x0
    ub[4] = Inf
    ub[5:6] = xf
    for i ∈ 1:3(N+1)
        if i > 2(N+1)
            lb[6 + 2*M + i] = - γ
            ub[6 + 2*M + i] = γ
        end
    end

    #X = zeros(2*(N+1))

    function con_fun(c, unk)
        offset = 6
        j = 1
        fx = Flow(f,variable=true)
        fvx = Flow(fv,autonomous=true,variable=false)
        buff = []
        X̃ = zeros(2*3)
        for (k,j) ∈ enumerate(J)
            if k == 1
                c[offset+1:offset+2] = [0.0,0.0] #x₀ = x̃₀
                push!(buff,offset+1,offset+2)
            else
                c[offset+2*k-1:offset+2*k] = unk[N+1+j+1:N+1+j+2] - X̃[5:6]
                push!(buff,offset+2*k-1:offset+2*k)
            end
            if j < J[end]
                for i ∈ 1:3
                    if i == 1
                        x_start = unk[N+1+j+1:N+1+j+2] # view(unk,N+1+j+1:N+1+j+2) #
                    else
                        x_start = X̃[2(i-1)-1:2(i-1)] # view(X̃,2(i-1)-1:2(i-1)) #
                    end
                    # tspan = (unk[j+i],unk[j+i+1])
                    # prob = ODEProblem(ff,x_start,tspan,unk[N+1 + 2*M + (i+j == N+1 ? i+j-1 : i + j)])
                    # sol = solve(prob,Tsit5())
                    # println(sol.u[end])
                    X̃[2i-1:2i] = fx(unk[j+i],x_start,unk[N+1 + 2*M + (i+j == N+1 ? i+j-1 : i + j)],unk[j+i+1]) 
                    #X̃[2i-1:2i] = fvx(unk[j+i],[x_start; unk[N+1 + 2*M + (i+j == N+1 ? i+j-1 : i + j)]],unk[j+i+1])[1:2]
                end
            end
        end
        offset = 6 + 2*M
        for i ∈ 1:2(N+1) # x
            c[i+offset] = 0
        end
        offset = 6 + 2*M
        for i ∈ 1:N+1 # u
            if i == N+1
                c[i+offset+2(N+1)] = unk[(N+1)+2*M+i-1]#U[end]
            else
                c[i+offset+2(N+1)] = unk[(N+1)+2*M+i]#U[i]
            end
        end
        c[1] = unk[1]
        c[2:3] = view(unk, (N+1)+1:(N+1)+2)
        c[4] = unk[N+1]
        c[5:6] = view(unk, (N+1)+2*M-1:(N+1)+2M)
        println(buff,c)
        return c
    end

    cx = similar(lb)
    con_fun(cx, unk0)
    adnlp = ADNLPModel!(obj_fun, unk0, l_var, u_var, con_fun, lb, ub)
    λa = [Float64(i) for i = 1:adnlp.meta.ncon]

    return adnlp, λa
end

adnlp, λa = solve_example()

solver_ad = IpoptSolver(adnlp)
ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol = 1e-12, mu_strategy="adaptive", sb="yes", print_level = 0)