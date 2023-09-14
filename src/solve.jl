"""
$(TYPEDSIGNATURES)

Solve a SimpleProblem with discretisation and return the solution : 
```(tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], var, obj```
"""
function solve(ocp::OptimalControlModel, N::Int64=9, M::Int64=(N + 1) ÷ 3; tol::Float64=1e-9, fixed_time_step::Bool=true)
  
  # check grids dimensions
  _check_grids(N,M)

  I = [i for i ∈ 0:N]
  J = 1:3:size(I, 1)

  initial_time_index = 1
  final_time_index = fixed_time_step ? 2 : N+1 
  time_size = final_time_index - initial_time_index + 1

  M = size(J, 1)

  (ξl, ξ, ξu), (ηl, η, ηu), (ψl, ψ, ψu), (ϕl, ϕ, ϕu), (θl, θ, θu), (ul, uind, uu), (xl, xind, xu), (vl, vind, vu) = nlp_constraints(ocp)
  
  parameter_dimension = ocp.control_dimension
  variable_dimension = is_variable_independent(ocp) ? 0 : ocp.variable_dimension
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
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.control_dimension,ocp.control_dimension+variable_dimension)] : zeros(eltype(arg),1)
    ẋ[:] = ocp.dynamics(t,x,arg[rg(1,ocp.control_dimension)],arg)
  end

  function dyna(x,arg,t)
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.control_dimension,ocp.control_dimension+variable_dimension)] : zeros(eltype(arg),1)
    return ocp.dynamics(t,x,arg[rg(1,ocp.control_dimension)],arg)
  end

  function f⁰(val,arg,t)
    x = arg[rg(1,ocp.state_dimension)]
    u = arg[rg(ocp.state_dimension+1,ocp.state_dimension+ocp.control_dimension)] 
    v = is_variable_dependent(ocp) ? arg[rg(1+ocp.state_dimension+ocp.control_dimension,ocp.state_dimension+ocp.control_dimension+variable_dimension)] : arg#zeros(eltype(arg),1)
    return ocp.lagrange(t,x,u,v)
  end

  # unknowns are tᵢ, xⱼ, βᵢ, λ, F
  unk_dim = final_time_index-initial_time_index+1 + ocp.state_dimension*M + parameter_dimension*N + variable_dimension + 1
  # constraints are ti+1 - ti, xj-̃xj, b, d(ₓ,ᵤ), c(ₓ,ᵤ,ₓᵤ), c(ᵥ), F
  con_dim = final_time_index-initial_time_index + ocp.state_dimension*M + bc_dimension + size(xind,1)*(N+1) + size(uind,1)*(N+1) + nlc_dimension*(N+1) + vc_dimension + 1 

  unk0 = _init_vector(N, 0.5, unk_dim)

  function obj_fun(unk)
    return(unk[end])# + norm(view(unk,1:N+1).-get_times_uniform(N,unk[1],unk[N+1])))
  end

  # l_var, u_var
  l_var = fill(-Inf,unk_dim)
  u_var = fill(Inf,unk_dim)
  for i ∈ initial_time_index:final_time_index # tᵢ constraints
      l_var[i] = 0.0
      u_var[i] = Inf
  end
  offset_unk = final_time_index - initial_time_index + 1

  # initial and final time if fixed
  if !CTBase.__is_initial_time_free(ocp)
    l_var[initial_time_index] = ocp.initial_time
    u_var[initial_time_index] = ocp.initial_time
  end
  if !CTBase.__is_final_time_free(ocp)
    l_var[final_time_index] = ocp.final_time
    u_var[final_time_index] = ocp.final_time
  end

  # box variable of ocp (!=unk) constraint
  offset_unk = offset_unk + ocp.state_dimension*M + parameter_dimension*N
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
  offset_con = fixed_time_step ? 1 : N 
  lb[1:offset_con] .= 10e-3
  ub[1:offset_con] .= Inf

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

    times = get_times(unk,N,fixed_time_step) # to change

    k = 1
    γ = variable_dimension > 0 ? unk[rg(time_size + M*ocp.state_dimension + N*parameter_dimension + 1,time_size + M*ocp.state_dimension + N*parameter_dimension + variable_dimension)] : Vector{eltype(unk)}()
    β = get_parameters(unk, ocp.state_dimension, parameter_dimension, N, M, time_size)
    lagrange_cost = 0
    xᵣ = ocp.state_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ocp.state_dimension)
    uᵢ = ocp.control_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ocp.control_dimension)
    
    if CTBase.__is_initial_time_free(ocp)
      variable_dimension > 1 ? γ[ocp.initial_time] = unk[initial_time_index] : γ = unk[initial_time_index] 
    end
    if CTBase.__is_final_time_free(ocp)
      variable_dimension > 1 ? γ[ocp.final_time] = unk[final_time_index] : γ = unk[final_time_index]
    end

    for i ∈ 1:N+1
      if i ∈ J
        xᵢ = unk[rg(time_size + (k-1)*ocp.state_dimension + 1, time_size + k*ocp.state_dimension)]
        if i == 1
          c[time_size-1+1:time_size-1+ocp.state_dimension] = zeros(eltype(unk),ocp.state_dimension)
        else
          c[rg(time_size-1 + (k-1)*ocp.state_dimension + 1, time_size-1 + k*ocp.state_dimension)] = xᵢ-xᵣ
        end
        k += 1
      else
        xᵢ = xᵣ
      end

      tᵢ = times[i]

      if i < N+1
        tᵢ₊₁ = times[i+1]
        if fixed_time_step
          c[1] = unk[final_time_index]-unk[initial_time_index]
        else
          c[i] = tᵢ₊₁ - tᵢ
        end
        tspan = (tᵢ,tᵢ₊₁)
        uᵢ = u_fun(β,i)
        arg = [uᵢ;γ]
        if ocp.state_dimension > 1
          p_dynamics = ODEProblem(dyna!, xᵢ, tspan, arg)
        else
          p_dynamics = ODEProblem(dyna, xᵢ, tspan, arg)
        end
        sol_dynamics = OrdinaryDiffEq.solve(p_dynamics, Tsit5(), abstol=1e-9)
        xᵣ = sol_dynamics(tᵢ₊₁)


        if !isnothing(ocp.lagrange)
          p_lagrange = ODEProblem(f⁰, 0.0, tspan, [xᵢ;uᵢ;γ])
          sol_lagrange = OrdinaryDiffEq.solve(p_lagrange, Tsit5(), abstol=1e-9)
          lagrange_cost += sol_lagrange(tᵢ₊₁)
        end


      end

      offset = time_size-1 + M*ocp.state_dimension + bc_dimension
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

    offset = time_size-1 + M*ocp.state_dimension
    if bc_dimension > 0
      c[offset+1:offset+bc_dimension] = ϕ(unk[rg(time_size + 1, time_size + ocp.state_dimension)], unk[rg(time_size + (M-1)*ocp.state_dimension + 1, time_size + M*ocp.state_dimension)], γ)
    end

    offset = offset + bc_dimension + (xind_dim + uind_dim + nlc_dimension)*(N+1)  
    if vc_dimension > 0
      c[offset+1:offset+vc_dimension] = θ(γ)
    end

    mayer_cost = 0
    if !isnothing(ocp.mayer)
      mayer_cost = ocp.mayer(unk[rg(time_size + 1, time_size + ocp.state_dimension)], unk[rg(time_size + (M-1)*ocp.state_dimension + 1, time_size + M*ocp.state_dimension)], γ)
    end

    c[end] = mayer_cost+lagrange_cost - unk[end]

    return c

  end

  cx = similar(lb)
  con_fun(cx, unk0)
  adnlp = ADNLPModel!(obj_fun, unk0, l_var, u_var, con_fun, lb, ub)

  solver_ad = IpoptSolver(adnlp)
  ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, adnlp, tol=1e-9, mu_strategy="adaptive", sb="yes", print_level=0)

  return(Solution(ipopt_solution_ad.solution,ocp,parameter_dimension,N,M,time_size))

end
