"""
$(TYPEDSIGNATURES)

Solve a SimpleProblem with discretisation and return the solution : 
```(tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], obj, var```
"""
function solve(prob::SimpleProblem, N::Int64=9, M::Int64=(N + 1) ÷ 3, tol::Float64=1e-9) #problem without variable λ and only boxconstraint on control and state

  TYPE::DataType = Float64

  # checking
  _check_grids(N, M)
  #M = M + 1
  function u(β, i)
    if i * prob.control_dim ≤ size(β, 1)
      val = view(β,prob.control_dim*(i-1)+1:prob.control_dim*i)
    else
      val = view(β,prob.control_dim*(i-2)+1:prob.control_dim*(i-1))
    end
    return size(prob.control_dim, 1) == 1 ? val[1] : val
  end

  function u!(val,β, i)
    if i * prob.control_dim > size(β, 1)
      i = i-1
    end
    for j in eachindex(val)
      val[j] = β[prob.control_dim*(i-1)+j] 
    end
  end

  I = [i for i ∈ 0:N]
  J = I[1:3:size(I, 1)]

  M = size(J, 1)

  # unknowns are (tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], obj, var  := unk
  function obj_fun(unk::Array{T}) where {T<:Real}
    # println("before get_times")
    # println("obj comp")
    #times = get_times(unk,N)
    #obj_m = unk[N+1+prob.state_dim*M+prob.control_dim*N+1] #+ norm(times.-get_times_uniform(prob,N,times[1],times[end]))
    #println(obj_m)
    return unk[N+1+prob.state_dim*M+prob.control_dim*N+1] #+ norm(times.-get_times_uniform(N,unk[1],unk[N+1])) #obj_m # obj + tᵢ constraints #unk[N+1+prob.state_dim*M+prob.control_dim*N+1] + norm(get_times(unk,N).-get_times_uniform(prob,N)) # obj + tᵢ constraints
  end

  dim_unk = N + 1 + prob.state_dim * M + prob.control_dim * N + 1 + 0 # tᵢ + xⱼ + βᵢ + F + λ

  unk0 = _init_vector(N, 0.5, dim_unk) #ones(Real,dim_unk)#

  l_var = fill(-Inf,dim_unk)
  u_var = fill(Inf,dim_unk)
  for i ∈ 1:dim_unk # variable bounds
    if i ≤ N + 1 # tᵢ constraints
      l_var[i] = 0.0
      u_var[i] = Inf
    end
  end
  l_var[end] = 0.0
  u_var[end] = Inf
  
  ff = ODEFunction{true}(prob.f!)

  # bf = BoundaryFreedom(prob)
  # if(!bf.x0)
  #   l_var[N+1+1:N+1+prob.state_dim] = prob.bmin[2:1+prob.state_dim]#[-1, 0]
  #   u_var[N+1+1:N+1+prob.state_dim] = prob.bmax[2:1+prob.state_dim]#[-1, 0]
  # end
  # if(!bf.xf)
  #   l_var[N+1+prob.state_dim*(M-1)+1:N+1+prob.state_dim*M] = prob.bmin[2+prob.state_dim+1:2+2*prob.state_dim]#[0, 0]
  #   u_var[N+1+prob.state_dim*(M-1)+1:N+1+prob.state_dim*M] = prob.bmin[2+prob.state_dim+1:2+2*prob.state_dim]#[0, 0]
  # end

  #println(l_var)
  #println(u_var)
  # the vector of constraint is :  t0, x0, tf, xf, xj-x̃j, x̃i, ui, tᵢ₊₁-tᵢ,F-obj
  dim_con = 2 * (1 + prob.state_dim) + M * prob.state_dim + (N + 1) * (prob.state_dim + prob.control_dim) + N + 1

  lb = zeros(dim_con)
  ub = zeros(dim_con)
  lb[1:2*(1+prob.state_dim)] = prob.bmin
  ub[1:2*(1+prob.state_dim)] = prob.bmax
  offset = 2 * (1 + prob.state_dim) + M * prob.state_dim
  for i ∈ 1:(prob.state_dim+prob.control_dim)*(N+1) # box constraint
    if i > prob.state_dim * (N + 1) && (i - 1) % prob.control_dim == 0# on u
      lb[range(offset + i, offset + i + prob.control_dim - 1)] = prob.dmin[range(prob.state_dim + 1, prob.state_dim + prob.control_dim)]
      ub[range(offset + i, offset + i + prob.control_dim - 1)] = prob.dmax[range(prob.state_dim + 1, prob.state_dim + prob.control_dim)]
    else # on x
      if (i - 1) % prob.state_dim == 0
        lb[range(offset + i, offset + i + prob.state_dim - 1)] = prob.dmin[range(1, prob.state_dim)]
        ub[range(offset + i, offset + i + prob.state_dim - 1)] = prob.dmax[range(1, prob.state_dim)]
      end
    end
  end
  # lb[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
  # ub[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
  # lb[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
  # ub[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
  # lb[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
  # ub[offset+1:offset+prob.state_dim] = [-1.0, 0.0]
  # lb[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
  # ub[offset+prob.state_dim*(M-1)+1:offset+prob.state_dim*(M-1)+prob.state_dim] = [0.0,0.0]
  lb[2*(1+prob.state_dim)+M*prob.state_dim+(N+1)*(prob.state_dim+prob.control_dim)+1:2*(1+prob.state_dim)+M*prob.state_dim+(N+1)*(prob.state_dim+prob.control_dim)+N] .= 1e-2
  ub[2*(1+prob.state_dim)+M*prob.state_dim+(N+1)*(prob.state_dim+prob.control_dim)+1:2*(1+prob.state_dim)+M*prob.state_dim+(N+1)*(prob.state_dim+prob.control_dim)+N] .= Inf
  lb[end] = 0
  ub[end] = 0

  function con_fun(c, unk)
    @assert size(lb, 1) == size(ub, 1) == size(c, 1)
    offset = 2 * (1 + prob.state_dim)
    j = 1
    #fx = Flow(prob.f, autonomous=false, variable=true)
    #buff = []
    X̃ = SizedVector{prob.state_dim * 3}(zeros(eltype(unk), prob.state_dim * 3)) # 3 is the number of state between two j
    β = SizedVector{prob.control_dim*N*1}(get_parameters(unk, prob.state_dim, prob.control_dim, N, M))#view(unk,N+1+prob.state_dim*M+1:N+1+prob.state_dim*M+N)
    l = SizedVector{prob.control_dim}(zeros(eltype(unk),prob.control_dim))
    zero_state = SizedVector{prob.state_dim}(zeros(eltype(unk),prob.state_dim))
    for (k, j) ∈ enumerate(J)
      if k == 1 # first index case
        c[offset+1:offset+prob.state_dim] = zero_state #x₀ = x̃₀  #unk[N+1+(k-1)*prob.state_dim+1:N+1+k*prob.state_dim]-[-1.0, 0.0]
        #push!(buff, offset+1:offset+prob.state_dim)
      else
        c[offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2)] = view(unk,N+1+(k-1)*prob.state_dim+1:N+1+k*prob.state_dim) - view(X̃,2*prob.state_dim+1:3*prob.state_dim)
        #push!(buff, offset+prob.state_dim*k-1:offset+prob.state_dim*k+(prob.state_dim-2))
      end     # to change when last index
      if j < J[end]
        for i ∈ 1:3
          #x_start = zeros(eltype(unk), prob.state_dim)
          if i == 1
            x_start = view(unk,N+1+prob.state_dim*(k-1)+1:N+1+prob.state_dim*k) # view(unk,N+1+j+1:N+1+j+2) #
          else
            x_start = view(X̃,prob.state_dim*(i-1)-1:prob.state_dim*(i-1)) # view(X̃,2(i-1)-1:2(i-1)) #
          end
          #println("x_start                 ",typeof(x_start))
          tspan = (unk[j+i], unk[j+i+1])
          # #println("tspan                   ",tspan)
          u!(l,β, i + j)
          #println("l                   ",typeof(l))
          p = ODEProblem(ff, x_start, tspan, l[rg(1,prob.control_dim)])
          sol = OrdinaryDiffEq.solve(p, Tsit5(), abstol=1e-9)
          X̃[prob.state_dim*(i-1)+1:prob.state_dim*i] = sol(unk[j+i+1])
          #println(sol(unk[j+i+1]))
          #X̃[prob.state_dim*(i-1)+1:prob.state_dim*i] = fx(unk[j+i], x_start, unk[j+i+1], u(β,i+j)) #   fx(0.0,x_start, 0.0, 1.0) #
        end
        c[offset+M*prob.state_dim+(k-1)*(3+1)*prob.state_dim+1:offset+M*prob.state_dim+k*(3+1)*prob.state_dim] = [unk[N+1+prob.state_dim*(k-1)+1:N+1+prob.state_dim*k]; X̃] # x̃ᵢ
      end
    end
    offset = 2 * (1 + prob.state_dim) + M * prob.state_dim
    # for i ∈ 1:prob.state_dim:prob.state_dim*(N+1) # x̃ᵢ 
    #     c[offset+i:offset+i+prob.state_dim-1] = zeros(prob.state_dim)
    #     push!(buff,offset+i:offset+i+prob.state_dim-1)
    # end
    lagrange_cost = 0.0
    for i ∈ 1:prob.control_dim:prob.control_dim*(N+1) # u
      if i == N + 1
        u!(view(c,offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1),β,i-1)#c[rg(offset+prob.state_dim*(N+1)+i,offset+prob.state_dim*(N+1)+i+prob.control_dim-1)] = u(β, i - 1)#unk[(N+1)+prob.state_dim*M+i-1]#U[end]
        #push!(buff, offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1)
      else
        u!(view(c,offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1),β,i)#c[rg(offset+prob.state_dim*(N+1)+i,offset+prob.state_dim*(N+1)+i+prob.control_dim-1)] = u(β, i)#unk[(N+1)+prob.state_dim*M+i]#U[i]
        function f⁰!(val,t)
          #ind = time_to_index(get_times(unk,N),t)

          val[1] = prob.f⁰(t,[0.0,0.0],0.0)#prob.f⁰(t,c[offset+(ind-1)*prob.state_dim+1:offset+ind*prob.state_dim],u(β, ind))
          return res
          #return t
        end
        f⁰x = Flow((t, _) -> prob.f⁰(t,[0.0,0.0],0.0), autonomous=false, variable=false)
        lagrange_cost += f⁰x((unk[i],unk[i+1]), zero(eltype(unk))).u[end]
        # #push!(buff, offset+prob.state_dim*(N+1)+i:offset+prob.state_dim*(N+1)+i+prob.control_dim-1)
      end
    end
    for i ∈ 1:N
      c[2*(1+prob.state_dim)+M*prob.state_dim+(N+1)*(prob.state_dim+prob.control_dim)+i] = unk[i+1] - unk[i]
    end
    c[1] = unk[1]
    c[2] = unk[N+1]
    c[3:4] = view(unk, (N+1)+1:(N+1)+prob.state_dim)
    c[5:6] = view(unk, (N+1)+2*M-(prob.state_dim-1):(N+1)+2M)
    mayer_cost = prob.g(unk[1], view(unk,N+2:N+3), unk[N+1], view(unk,N+1+(M-1)*prob.state_dim+1:N+1+M*prob.state_dim))
    c[end] = mayer_cost+lagrange_cost - unk[end]

    #push!(buff, 1:6)
    #println(buff)
    #println(c)
    return c
  end
  #replace!() TOSEE  
  #map!() TOSEE
  cx = similar(lb)
  con_fun(cx, unk0)
  adnlp = ADNLPModel!(obj_fun, unk0, l_var, u_var, con_fun, lb, ub)

  solver_ad = IpoptSolver(adnlp)
  ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol=1e-9, mu_strategy="adaptive", sb="yes", print_level=5)

  return (Solution(ipopt_solution_ad.solution,prob,N,M))

end


function solve(ocp::OptimalControlModel, N::Int64=9, M::Int64=(N + 1) ÷ 3, tol::Float64=1e-9, fixed_time_step::Bool=true)
  
  # check grids dimensions
  _check_grids(N,M)

  I = [i for i ∈ 0:N]
  J = 1:3:size(I, 1)

  M = size(J, 1)

  (ξl, ξ, ξu), (ηl, η, ηu), (ψl, ψ, ψu), (ϕl, ϕ, ϕu), (θl, θ, θu), (ul, uind, uu), (xl, xind, xu), (vl, vind, vu) = nlp_constraints(ocp)
  
  parameter_dimension = ocp.control_dimension
  nlc_dimension = size(ηl,1) + size(ξl,1) + size(ψl,1)
  vc_dimension = size(θl,1) 
  bc_dimension = size(ϕl,1)

  function u!(val,β, i)
    if i * ocp.control_dimension > size(β, 1)
      i = i-1
    end
    for j in eachindex(val)
      val[j] = β[ocp.control_dimension*(i-1)+j] 
    end
  end

  function u_fun(β, i)
    if i * ocp.control_dimension ≤ size(β, 1)
      val = view(β,ocp.control_dimension*(i-1)+1:ocp.control_dimension*i)
    else
      val = view(β,ocp.control_dimension*(i-2)+1:ocp.control_dimension*(i-1))
    end
    return ocp.control_dimension == 1 ? val[1] : val
  end

  function dyna!(ẋ,x,arg,t)
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.control_dimension,ocp.control_dimension+ocp.variable_dimension)] : zero(eltype(x))
    ẋ[:] = ocp.dynamics(t,x,arg[rg(1,ocp.control_dimension)],v)
  end

  function dyna(x,arg,t)
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.control_dimension,ocp.control_dimension+ocp.variable_dimension)] : zero(eltype(x))
    return ocp.dynamics(t,x,arg[rg(1,ocp.control_dimension)],v)
  end

  function f⁰(val,arg,t)
    x = arg[rg(1,ocp.state_dimension)]
    u = arg[rg(ocp.state_dimension+1,ocp.state_dimension+ocp.control_dimension)] 
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.state_dimension+ocp.control_dimension,ocp.state_dimension+ocp.control_dimension+ocp.variable_dimension)] : zero(eltype(x))
    return ocp.lagrange(t,x,u,v)
  end

  # unknowns are tᵢ, xⱼ, βᵢ, λ, F
  unk_dim = N+1 + ocp.state_dimension*M + parameter_dimension*N + ocp.variable_dimension + 1
  # constraints are ti+1 - ti, xj-̃xj, b, d(ₓ,ᵤ), c(ₓ,ᵤ,ₓᵤ), c(ᵥ), F
  con_dim = N + ocp.state_dimension*M + bc_dimension + size(xind,1)*(N+1) + size(uind,1)*(N+1) + nlc_dimension*(N+1) + vc_dimension + 1 

  unk0 = _init_vector(N, 0.5, unk_dim)

  function obj_fun(unk::Array{T}) where {T<:Real}
    return(unk[end])# + norm(view(unk,1:N+1).-get_times_uniform(N,unk[1],unk[N+1])))
  end

  # l_var, u_var
  l_var = fill(-Inf,unk_dim)
  u_var = fill(Inf,unk_dim)
  for i ∈ 1:N+1 # tᵢ constraints
      l_var[i] = 0.0
      u_var[i] = Inf
  end
  
  # initial and final time if fixed
  if !CTBase.__is_initial_time_free(ocp)
    l_var[1] = ocp.initial_time
    u_var[1] = ocp.initial_time
  end
  if !CTBase.__is_final_time_free(ocp)
    l_var[N+1] = ocp.final_time
    u_var[N+1] = ocp.final_time
  end

  # box variable of ocp (!=unk) constraint
  offset_unk = N+1 + ocp.state_dimension*M + parameter_dimension*N
  if size(vind,1) > 0
    for (i,v) in enumerate(vind)
      l_var[offset_unk+v] = vl[i]
      u_var[offset_unk+v] = vu[i]
    end
  end

  # cost constraint (maybe to remove)
  l_var[end] = 0.0
  u_var[end] = Inf
  
  # lb, ub
  lb = zeros(con_dim)
  ub = zeros(con_dim)

  # ti constraints
  lb[1:N] .= 10e-3
  ub[1:N] .= Inf
  offset_con = N 

  # boundary constraints
  offset_con = offset_con + ocp.state_dimension*M
  if bc_dimension > 0
    lb[offset_con+1:offset_con+bc_dimension] = ϕl
    ub[offset_con+1:offset_con+bc_dimension] = ϕu
  end

  # box state and control 
  offset_con = offset_con + bc_dimension
  xind_dim = size(xind,1)
  if xind_dim > 0
    for i ∈ 1:N+1
      lb[offset_con + (i-1)*xind_dim + 1 : offset_con + i*xind_dim] = xl
      ub[offset_con + (i-1)*xind_dim + 1 : offset_con + i*xind_dim] = xu
    end
  end
  offset_con =  offset_con + xind_dim*(N+1)
  uind_dim = size(uind,1)
  if uind_dim > 0
    for i ∈ 1:N+1
      lb[offset_con + (i-1)*uind_dim + 1 : offset_con + i*uind_dim] = ul
      ub[offset_con + (i-1)*uind_dim + 1 : offset_con + i*uind_dim] = uu
    end
  end

  # non linear constraints
  offset_con =  offset_con + uind_dim*(N+1)
  if nlc_dimension > 0
    nlc(t,x,u,v) = [η(t,x,v);ξ(t,u,v);ψ(t,x,u,v)]
    nlcl = [ηl;ξl;ψl]
    nlcu = [ηu;ξu;ψu]
    for i ∈ 1:N+1
      lb[offset_con + (i-1)*nlc_dimension + 1 : offset_con + i*nlc_dimension] = nlcl
      ub[offset_con + (i-1)*nlc_dimension + 1 : offset_con + i*nlc_dimension] = nlcu
    end
  end
  
  offset_con = offset_con + (N+1)*nlc_dimension
  if vc_dimension > 1
    lb[offset_con + 1 : offset_con + vc_dimension] = θl
    ub[offset_con + 1 : offset_con + vc_dimension] = θu
  end


  function con_fun(c,unk)
    @assert size(lb, 1) == size(ub, 1) == size(c, 1)

    offset = 0
    k = 1
    γ = ocp.variable_dimension > 0 ? unk[rg(N+1 + M*ocp.state_dimension + N*parameter_dimension + 1,N+1 + M*ocp.state_dimension + N*parameter_dimension + ocp.variable_dimension)] : []
    β = get_parameters(unk, ocp.state_dimension, parameter_dimension, N, M)
    lagrange_cost = 0
    xᵣ = ocp.state_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ocp.state_dimension)
    uᵢ = ocp.control_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ocp.control_dimension)
    
    if fixed_time_step # change time between 0 and 1 (* tf-t0)
      unk[1:N+1] = get_times_uniform(N,unk[1],unk[N+1])
    end
    
    if CTBase.__is_initial_time_free(ocp)
      ocp.variable_dimension > 1 ? γ[ocp.initial_time] = unk[1] : γ = unk[1] 
    end
    if CTBase.__is_final_time_free(ocp)
      ocp.variable_dimension > 1 ? γ[ocp.final_time] = unk[N+1] : γ = unk[N+1]
    end

    for i ∈ 1:N+1
      if i ∈ J
        xᵢ = unk[rg(N+1 + (k-1)*ocp.state_dimension + 1, N+1 + k*ocp.state_dimension)]
        if i == 1
          c[N+1:N+ocp.state_dimension] = zeros(eltype(unk),ocp.state_dimension)
        else
          c[rg(N + (k-1)*ocp.state_dimension + 1, N + k*ocp.state_dimension)] = xᵢ-xᵣ
        end
        k += 1
      else
        xᵢ = xᵣ
      end

      if i < N+1
        tᵢ = unk[i]
        tᵢ₊₁ = unk[i+1]
        c[i] = tᵢ₊₁ - tᵢ
        tspan = (tᵢ,tᵢ₊₁)
        uᵢ = u_fun(β,i)
        arg = [uᵢ;γ]
        if ocp.state_dimension > 1
          p_dynamics = ODEProblem(dyna!, xᵢ, tspan, arg)
        else
          # println(dyna)
          # println(xᵢ)
          # println(tspan)
          # println(arg)
          p_dynamics = ODEProblem(dyna, xᵢ, tspan, arg)
        end
        sol_dynamics = OrdinaryDiffEq.solve(p_dynamics, Tsit5(), abstol=1e-9)
        xᵣ = sol_dynamics(tᵢ₊₁)


        if !isnothing(ocp.lagrange)
          p_lagrange = ODEProblem(f⁰, lagrange_cost, tspan, [xᵢ;uᵢ;γ])
          sol_lagrange = OrdinaryDiffEq.solve(p_lagrange, Tsit5(), abstol=1e-9)
          lagrange_cost += sol_lagrange(tᵢ₊₁)
        end


      end

      offset = N + M*ocp.state_dimension + bc_dimension
      if xind_dim > 0
        c[offset + (i-1)*xind_dim + 1 : offset + i*xind_dim] = [xᵢ[ind] for ind ∈ xind]
      end
      
      offset = offset + (N+1)*xind_dim
      if uind_dim > 0
        c[offset + (i-1)*uind_dim + 1 : offset + i*uind_dim] = [uᵢ[ind] for ind ∈ uind]
      end
      
      offset = offset + (N+1)*uind_dim
      if nlc_dimension > 0
        c[offset + (i-1)*nlc_dimension + 1 : offset + i*nlc_dimension] = nlc(tᵢ,xᵢ,uᵢ,γ)
      end

    end

    offset = N + M*ocp.state_dimension
    if bc_dimension > 0
      c[offset+1:offset+bc_dimension] = ϕ(unk[rg(N+1 + 1, N+1 + ocp.state_dimension)], unk[rg(N+1 + (M-1)*ocp.state_dimension + 1, N+1 + M*ocp.state_dimension)], γ)
    end

    offset = offset + bc_dimension + (xind_dim + uind_dim + nlc_dimension)*(N+1)  
    if vc_dimension > 0
      c[offset+1:offset+vc_dimension] = θ(γ)
    end

    mayer_cost = 0
    if !isnothing(ocp.mayer)
      mayer_cost = ocp.mayer(unk[rg(N+1 + 1, N+1 + ocp.state_dimension)], unk[rg(N+1 + (M-1)*ocp.state_dimension + 1, N+1 + M*ocp.state_dimension)], γ)
    end

    c[end] = mayer_cost+lagrange_cost - unk[end]

    return c

  end

  println(l_var)
  println(u_var)

  println(lb)
  println(ub)

  cx = similar(lb)
  con_fun(cx, unk0)
  adnlp = ADNLPModel!(obj_fun, unk0, l_var, u_var, con_fun, lb, ub)

  solver_ad = IpoptSolver(adnlp)
  ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol=1e-9, mu_strategy="adaptive", sb="yes", print_level=5)

  return(Solution(ipopt_solution_ad.solution,ocp,parameter_dimension,N,M))

end
