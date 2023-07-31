function solve(prob::SimpleProblem, N::Int64=9, M::Int64=(N+1)÷3) #problem without variable λ and only boxconstraint on control and state
    
    TYPE::DataType = Float64

    # checking
    _check_grids(N,M)
    M = M+1
    function u(β,i) 
        if i*prob.control_dim ≤ size(β,1)
            val = prob.l(β[prob.control_dim*(i-1)+1:prob.control_dim*i])
        else
            val = prob.l(β[prob.control_dim*(i-2)+1:prob.control_dim*(i-1)])
        end
        return size(val,1) == 1 ? val[1] : val  
    end

    # make the dynamics in-place TODO
    function f!(dx::Array{T},x::Array{Y},u,t) where {T<:Real, Y<:Real}
        dx[1] = x[2]
        dx[2] = u
    end

    I = [i for i ∈ 0:N]
    J = I[1:3:size(I,1)]
    
    M = size(J,1)

    # unknowns are (tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], obj, var)  := unk
    function obj_fun(unk::Array{T})  where T<:Real
        return unk[N+1] + norm(get_times(unk,N).-get_times_uniform(prob,N)) # obj + tᵢ constraints #unk[N+1+prob.state_dim*M+prob.control_dim*N+1] + norm(get_times(unk,N).-get_times_uniform(prob,N)) # obj + tᵢ constraints
    end

    dim_unk =  N+1 + prob.state_dim*M + prob.control_dim*N + 1 + 0 # tᵢ + xⱼ + βᵢ + F + λ
    
    unk0 = _init_vector(N,0.5,dim_unk) #ones(Real,dim_unk)#

    l_var = Vector{TYPE}(-Inf*ones(dim_unk))
    u_var = Vector{TYPE}(Inf*ones(dim_unk))
    for i ∈ 1:dim_unk # variable bounds
        if i ≤ N+1 # tᵢ constraints
            l_var[i] = 0.0
            u_var[i] = Inf
        end
    end
    l_var[end] = 0.0
    u_var[end] = Inf
    l_var[N+1+1:N+1+prob.state_dim] = [-1,0]
    u_var[N+1+1:N+1+prob.state_dim] = [-1,0]
    l_var[N+1+prob.state_dim*(M-1)+1:N+1+prob.state_dim*M] = [0,0]
    u_var[N+1+prob.state_dim*(M-1)+1:N+1+prob.state_dim*M] = [0,0]

    println(l_var)
    println(u_var)
    # the vector of constraint is :  t0, x0, tf, xf, xj-x̃j, x̃i, ui, λ
    dim_con = 2*(1 + prob.state_dim) + M*prob.state_dim + (N+1)*(prob.state_dim + prob.control_dim) + 1

    lb = zeros(TYPE, dim_con)
    ub = zeros(TYPE, dim_con)
    lb[1:2*(1 + prob.state_dim)] = prob.bmin
    ub[1:2*(1 + prob.state_dim)] = prob.bmax
    offset = 2*(1 + prob.state_dim) + M*prob.state_dim
    for i ∈ 1:(prob.state_dim + prob.control_dim)*(N+1) # box constraint
        if i > prob.state_dim*(N+1) && (i-1)%prob.control_dim==0# on u
            lb[range(offset + i,offset + i + prob.control_dim-1)] = prob.dmin[range(prob.state_dim + 1,prob.state_dim + prob.control_dim)]
            ub[range(offset + i,offset + i + prob.control_dim-1)] = prob.dmax[range(prob.state_dim + 1,prob.state_dim + prob.control_dim)]
        else # on x
            if (i-1)%prob.state_dim==0
                lb[range(offset+i,offset+i+prob.state_dim-1)] = prob.dmin[range(1,prob.state_dim)]
                ub[range(offset+i,offset+i+prob.state_dim-1)] = prob.dmax[range(1,prob.state_dim)]
            end
        end
    end
    # lb[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
    # ub[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
    # lb[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
    # ub[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
    lb[end] = 0
    ub[end] = Inf

    function con_fun(c::Array{T}, unk::Array{Y}) where {T<:Real, Y<:Real}
        @assert size(lb,1) == size(ub,1) == size(c,1)
        offset = 2*(1 + prob.state_dim)
        j = 1
        fx = Flow(prob.f, autonomous=false, variable=true)
        buff = []
        X̃ = zeros(prob.state_dim*3) # 3 is the number of state between to j
        β = get_parameters(unk,prob.state_dim,prob.control_dim,N,M)#view(unk,N+1+prob.state_dim*M+1:N+1+prob.state_dim*M+N)
        for (k,j) ∈ enumerate(J)
            if k == 1 # first index case
                c[offset+1:offset+prob.state_dim] = zeros(prob.state_dim) #x₀ = x̃₀ 
                push!(buff,offset+1:offset+prob.state_dim)
            else
                c[offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2)] = unk[N+1+j+1:N+1+j+2] - X̃[2*prob.state_dim+1:3*prob.state_dim]
                push!(buff,offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2))
            end     # to change when last index
            if j < J[end]
                for i ∈ 1:3
                    x_start = Array{eltype(unk)}(undef,prob.state_dim)
                    if i == 1
                        x_start = unk[N+1+j+1:N+1+j+prob.state_dim] # view(unk,N+1+j+1:N+1+j+2) #
                    else
                        x_start = X̃[prob.state_dim*(i-1)-1:prob.state_dim*(i-1)] # view(X̃,2(i-1)-1:2(i-1)) #
                    end
                    #println("x_start                 ",typeof(x_start))
                    tspan = (unk[j+i],unk[j+i+1])
                    #println("tspan                   ",tspan)
                    l = u(β,i+j)
                    #println("l                   ",typeof(l))
                    p = ODEProblem(f!,Array{Float64}(x_start),tspan,l)
                    sol = OrdinaryDiffEq.solve(p,Tsit5())
                    X̃[prob.state_dim*(i-1)+1:prob.state_dim*i] = Array{eltype(unk)}(sol(unk[j+i+1]))
                    println(sol(unk[j+i+1]))
                    #X̃[prob.state_dim*(i-1)+1:prob.state_dim*i] =fx(unk[j+i], x_start, unk[j+i+1], u(β,i+j)) #   fx(0.0,x_start, 0.0, 1.0) #
                end
            end
        end
        offset = 2*(1 + prob.state_dim) + M*prob.state_dim
        for i ∈ 1:prob.state_dim:prob.state_dim*(N+1) # x̃ᵢ 
            c[offset+i:offset+i+prob.state_dim-1] = zeros(prob.state_dim)
            push!(buff,offset+i:offset+i+prob.state_dim-1)
        end
        for i ∈ 1:prob.control_dim:prob.control_dim*(N+1) # u
            if i == N+1
                c[offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1] = [u(β, i-1)]#unk[(N+1)+prob.state_dim*M+i-1]#U[end]
                push!(buff,offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1)
            else
                c[offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1] = [u(β, i)]#unk[(N+1)+prob.state_dim*M+i]#U[i]
                push!(buff,offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1)
            end
        end
        c[1] = unk[1]
        c[2:3] = view(unk, (N+1)+1:(N+1)+prob.state_dim)
        c[4] = unk[N+1]
        c[5:6] = view(unk, (N+1)+2*M-(prob.state_dim-1):(N+1)+2M)
        c[end] = prob.g(unk[1],unk[N+2:N+3],unk[N+1],unk[N+1+M-1:N+1+M])

        push!(buff,1:6)
        #println(buff)
        #println(c)
        return c
    end

    cx = similar(lb)
    con_fun(cx, unk0)
    adnlp = ADNLPModel!(obj_fun, unk0, l_var, u_var, con_fun, lb, ub)
    λa = [Float64(i) for i = 1:adnlp.meta.ncon]

    x, fx, ngx, optimal, iter = steepest(adnlp)
    println("x = $x")
    println("fx = $fx")
    println("ngx = $ngx")
    println("optimal = $optimal")
    println("iter = $iter")
    # solver_ad = IpoptSolver(adnlp)
    # ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol = 1e-12, mu_strategy="adaptive", sb="yes", print_level = 0)
    
    return(true)

end


# ref : https://discourse.julialang.org/t/dual-forwarddiff-type-error-when-jump-tries-autodiff-user-defined-function/25639