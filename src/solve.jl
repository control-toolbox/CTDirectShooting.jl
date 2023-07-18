function solve(prob::SimpleProblem, N::Int64=9, M::Int64=N÷3) #problem without variable λ and only boxconstraint on control and state
    
    # checking
    _check_grids(N,M)

    # transform functions for scalar/vectors
    dyn(x,u) = size(u) == 1 ? prob.f(x,u[1]) : prob.f(x,u)

    boundary_freedom = BoundaryFreedom(prob)

    #u(β,i) = i*prob.control_dim ≤ size(β,1) ? prob.l(β[prob.control_dim*(i-1)+1:prob.control_dim*i]) : prob.l(β[prob.control_dim*(i-2)+1:prob.control_dim*(i-1)]) # picewise constant !! DIM 1 /other todo
    u(β,i) = i ≤ size(β,1) ? prob.l(β[i]) : prob.l(β[i-1]) # picewise constant !! DIM 1 /other todo


    I = [i for i ∈ 0:N]
    J = I[1:3:size(I,1)]
    
    # unknowns are (tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], obj, var)  := unk
    obj_fun(unk) = unk[N+1+M+N+1] + norm(get_times(unk).-get_times_uniform(prob)) # obj + tᵢ constraints
    
    dim_unk =  N+1 + 2*M + 1*N + 1 + 0 # tᵢ + xⱼ + βᵢ + F + λ
    
    unk0 = _init_vector(dim_unk,0.5)

    l_var = -Inf*ones(dim_unk)
    u_var = Inf*ones(dim_unk)
    for i ∈ 1:dim_unk # variable bounds
        if i ≤ N+1 # tᵢ constraints
            l_var[i] = 0
            u_var[i] = Inf
        end
    end

    # the vector of constraint is :  t0,x0,tf,xf,
    dim_con = 2*(1 + prob.state_dim) + M*prob.state_dim + (N+1)*(prob.state_dim + prob.control_dim) 

    lb = zeros(dim_con)
    ub = zeros(dim_con)
    lb[1:2*(1 + prob.state_dim)] = prob.bmin
    ub[1:2*(1 + prob.state_dim)] = prob.bmax
    offset = 2*(1 + prob.state_dim) + M*prob.state_dim
    for i ∈ 1:(prob.state_dim + prob.control_dim)*(N+1) # box constraint
        if i > prob.state_dim*(N+1) && i%prob.control_dim==1# on u
            lb[range(offset + i,offset + i + prob.control_dim-1)] = prob.dmin[range(prob.state_dim + 1,prob.state_dim + prob.control_dim)]
            ub[range(offset + i,offset + i + prob.control_dim-1)] = prob.dmax[range(prob.state_dim + 1,prob.state_dim + prob.control_dim)]
        else # on x
            if i%prob.state_dim==1
                lb[range(offset+i,offset+i+prob.state_dim-1)] = prob.dmin[range(1,prob.state_dim)]
                ub[range(offset+i,offset+i+prob.state_dim-1)] = prob.dmax[range(1,prob.state_dim)]
            end
        end
    end

    function con_fun(c, unk)
        offset = 2*(1 + prob.state_dim) + M*prob.state_dim
        j = 1
        fx = Flow(prob.f, variable=true)
        buff = []
        X̃ = zeros(prob.state_dim*3) # 3 is the number of state between to j
        for (k,j) ∈ enumerate(J)
            if k == 1
                c[offset+1:offset+prob.state_dim] = zeros(prob.state_dim) #x₀ = x̃₀
                push!(buff,offset+1:offset+prob.state_dim)
            else
                c[offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2)] = unk[N+1+j+1:N+1+j+2] - X̃[2*prob.state_dim+1:3*prob.state_dim]
                push!(buff,offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2))
            end
            if j < J[end]
                for i ∈ 1:3
                    if i == 1
                        x_start = unk[N+1+j+1:N+1+j+prob.state_dim] # view(unk,N+1+j+1:N+1+j+2) #
                    else
                        x_start = X̃[prob.state_dim*(i-1)-1:prob.state_dim*(i-1)] # view(X̃,2(i-1)-1:2(i-1)) #
                    end
                    X̃[prob.state_dim*(i-1)+1:prob.state_dim*i] = fx(unk[j+i],x_start, u(unk[N+1+prob.state_dim*M+1:end],i+j), unk[j+i+1]) # replace with view 
                end
            end
        end # todo
        for i ∈ 1:prob.state_dim*(N+1) # x /
            c[i+offset] = 0
        end
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

    solver_ad = IpoptSolver(adnlp)
    ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol = 1e-12, mu_strategy="adaptive", sb="yes", print_level = 0)
    
    return(ipopt_solution_ad)

end