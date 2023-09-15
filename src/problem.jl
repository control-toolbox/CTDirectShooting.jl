"""
$(TYPEDEF)

**Fields**

$(TYPEDFIELDS)

"""
mutable struct CTDirectShooting_data

  ocp::OptimalControlModel
  grid_size_fine::Integer
  grid_size_coarse::Integer
  
  unk_dim
  con_dim
  initial_time_index
  final_time_index
  time_size

  fixed_time_step

  parameter_dimension
  variable_dimension

  # nlp_constraints
  control_constraints
  state_constraints
  mixed_constraints
  boundary_conditions
  variable_constraints
  control_box
  state_box
  variable_box

  nlc
  nlc_dimension
  vc_dimension
  bc_dimension

  u_fun
  dyna!
  dyna
  f⁰

  init
  tol

  function CTDirectShooting_data(ocp,grid_size_fine,grid_size_coarse,fixed_time_step,init,tol)
    
    ctds = new()
    
    ctds.ocp = ocp
    ctds.grid_size_fine = grid_size_fine
    ctds.grid_size_coarse = grid_size_coarse

    ctds.initial_time_index = 1
    ctds.final_time_index = fixed_time_step ? 2 : grid_size_fine+1 
    ctds.time_size = ctds.final_time_index - ctds.initial_time_index + 1

    ctds.fixed_time_step = fixed_time_step

    ctds.parameter_dimension = ctds.ocp.control_dimension
    ctds.variable_dimension = is_variable_independent(ctds.ocp) ? 0 : ctds.ocp.variable_dimension
    ctds.unk_dim = ctds.final_time_index-ctds.initial_time_index+1 + ctds.ocp.state_dimension*(ctds.grid_size_coarse+1) + ctds.parameter_dimension*grid_size_fine + ctds.variable_dimension + 1
    
    ctds.control_constraints, ctds.state_constraints, ctds.mixed_constraints, ctds.boundary_conditions, ctds.variable_constraints, ctds.control_box, ctds.state_box, ctds.variable_box = nlp_constraints(ocp)

    ctds.nlc = (t,x,u,v) -> [(ctds.state_constraints[2])(t,x,v);(ctds.control_constraints[2])(t,u,v);(ctds.mixed_constraints[2])(t,x,u,v)]
    ctds.nlc_dimension = size(ctds.state_constraints[1],1) + size(ctds.control_constraints[1],1) + size(ctds.mixed_constraints[1],1)
    ctds.vc_dimension = size(ctds.variable_constraints[1],1) 
    ctds.bc_dimension = size(ctds.boundary_conditions[1],1)

    ctds.con_dim = ctds.final_time_index-ctds.initial_time_index + ctds.ocp.state_dimension*(ctds.grid_size_coarse+1) + ctds.bc_dimension + size(ctds.state_box[2],1)*(ctds.grid_size_fine+1) + size(ctds.control_box[2],1)*(ctds.grid_size_fine+1) + ctds.nlc_dimension*(ctds.grid_size_fine+1) + ctds.vc_dimension + 1 

    ctds.u_fun = (β,i) -> begin
      if i * ctds.ocp.control_dimension ≤ size(β, 1)
        val = view(β,ctds.ocp.control_dimension*(i-1)+1:ctds.ocp.control_dimension*i)
      else
        val = view(β,ctds.ocp.control_dimension*(i-2)+1:ctds.ocp.control_dimension*(i-1))
      end
      return ctds.ocp.control_dimension == 1 ? val[1] : val
    end

    ctds.dyna! = (ẋ,x,arg,t) -> begin
      v = is_variable_dependent(ctds.ocp) ? arg[rg(1+ctds.ocp.control_dimension,ctds.ocp.control_dimension+ctds.variable_dimension)] : arg
      ẋ[:] = ctds.ocp.dynamics(t,x,arg[rg(1,ctds.ocp.control_dimension)],v)
    end
  
    ctds.dyna = (x,arg,t) -> begin
      v = is_variable_dependent(ctds.ocp) ? arg[rg(1+ctds.ocp.control_dimension,ctds.ocp.control_dimension+ctds.variable_dimension)] : arg
      return ctds.ocp.dynamics(t,x,arg[rg(1,ctds.ocp.control_dimension)],v)
    end

    ctds.f⁰ = (val,arg,t) -> begin
      x = arg[rg(1,ctds.ocp.state_dimension)]
      u = arg[rg(ctds.ocp.state_dimension+1,ctds.ocp.state_dimension+ctds.ocp.control_dimension)] 
      v = is_variable_dependent(ctds.ocp) ? arg[rg(1+ctds.ocp.state_dimension+ctds.ocp.control_dimension,ctds.ocp.state_dimension+ctds.ocp.control_dimension+ctds.variable_dimension)] : arg
      return ctds.ocp.lagrange(t,x,u,v)
    end

    ctds.init = init
    ctds.tol = tol

    return ctds
  end

end

# compute l_var and u_var
function variable_bounds(ctds)

  vl, vind, vu = ctds.variable_box

  # l_var, u_var
  l_var = fill(-Inf,ctds.unk_dim)
  u_var = fill(Inf,ctds.unk_dim)
  for i ∈ ctds.initial_time_index:ctds.final_time_index # tᵢ constraints
    l_var[i] = 0.0
    u_var[i] = Inf
  end
  offset_unk = ctds.final_time_index - ctds.initial_time_index + 1

  # initial and final time if fixed
  if !CTBase.__is_initial_time_free(ctds.ocp)
    l_var[ctds.initial_time_index] = ctds.ocp.initial_time
    u_var[ctds.initial_time_index] = ctds.ocp.initial_time
  end
  if !CTBase.__is_final_time_free(ctds.ocp)
    l_var[ctds.final_time_index] = ctds.ocp.final_time
    u_var[ctds.final_time_index] = ctds.ocp.final_time
  end

  # box variable of ocp (!=unk) constraint
  offset_unk = offset_unk + ctds.ocp.state_dimension*(ctds.grid_size_coarse+1) + ctds.parameter_dimension*ctds.grid_size_fine
  if size(vind,1) > 0
    for (i,v) in enumerate(vind)
      l_var[offset_unk+v] = vl[i]
      u_var[offset_unk+v] = vu[i]
    end
  end

  # cost constraint (maybe to remove)
  l_var[end] = 0.0
  u_var[end] = Inf

  return l_var,u_var

end

function constraint_bounds(ctds)
  # lb, ub
  lb = zeros(ctds.con_dim)
  ub = zeros(ctds.con_dim)

  # ti constraints
  offset_con = ctds.fixed_time_step ? 1 : ctds.grid_size_fine 
  lb[1:offset_con] .= 10e-3
  ub[1:offset_con] .= Inf

  # boundary conditions
  offset_con = offset_con + ctds.ocp.state_dimension*(ctds.grid_size_coarse+1)
  ϕl, _, ϕu = ctds.boundary_conditions 
  if ctds.bc_dimension > 0
    lb[offset_con+1:offset_con+ctds.bc_dimension] = ϕl
    ub[offset_con+1:offset_con+ctds.bc_dimension] = ϕu
  end

  # box state and control 
  offset_con = offset_con + ctds.bc_dimension
  xl, xind, xu = ctds.state_box
  xind_dim = size(xind,1)
  if xind_dim > 0
    for i ∈ 1:ctds.grid_size_fine+1
      lb[offset_con + (i-1)*xind_dim + 1 : offset_con + i*xind_dim] = xl
      ub[offset_con + (i-1)*xind_dim + 1 : offset_con + i*xind_dim] = xu
    end
  end
  offset_con =  offset_con + xind_dim*(ctds.grid_size_fine+1)
  ul, uind, uu = ctds.control_box
  uind_dim = size(uind,1)
  if uind_dim > 0
    for i ∈ 1:ctds.grid_size_fine+1
      lb[offset_con + (i-1)*uind_dim + 1 : offset_con + i*uind_dim] = ul
      ub[offset_con + (i-1)*uind_dim + 1 : offset_con + i*uind_dim] = uu
    end
  end

  # non linear constraints
  offset_con =  offset_con + uind_dim*(ctds.grid_size_fine+1)
  ηl, _, ηu = ctds.state_constraints
  ξl, _, ξu = ctds.control_constraints
  ψl, _, ψu = ctds.mixed_constraints
  if ctds.nlc_dimension > 0
    #nlc(t,x,u,v) = [η(t,x,v);ξ(t,u,v);ψ(t,x,u,v)]
    nlcl = [ηl;ξl;ψl]
    nlcu = [ηu;ξu;ψu]
    for i ∈ 1:ctds.grid_size_fine+1
      lb[offset_con + (i-1)*ctds.nlc_dimension + 1 : offset_con + i*ctds.nlc_dimension] = nlcl
      ub[offset_con + (i-1)*ctds.nlc_dimension + 1 : offset_con + i*ctds.nlc_dimension] = nlcu
    end
  end
  
  offset_con = offset_con + (ctds.grid_size_fine+1)*ctds.nlc_dimension
  if ctds.vc_dimension > 1
    lb[offset_con + 1 : offset_con + ctds.vc_dimension] = θl
    ub[offset_con + 1 : offset_con + ctds.vc_dimension] = θu
  end

  return lb,ub
end

function obj_fun(unk,ctds)
  return(unk[end])
end

function con_fun(c,unk,ctds)

  times = get_times(unk,ctds)

  k = 1
  γ = ctds.variable_dimension > 0 ? get_variable(unk,ctds) : Vector{eltype(unk)}()
  β = get_parameters(unk, ctds)
  lagrange_cost = 0
  xᵣ = ctds.ocp.state_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ctds.ocp.state_dimension)
  uᵢ = ctds.ocp.control_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ctds.ocp.control_dimension)
  
  _, xind, _ = ctds.state_box
  _, uind, _ = ctds.control_box
  _, ϕ, _ = ctds.boundary_conditions
  _, θ, _ = ctds.variable_constraints
  xind_dim = size(xind,1)
  uind_dim = size(uind,1)

  if CTBase.__is_initial_time_free(ctds.ocp)
    ctds.variable_dimension > 1 ? γ[ctds.ocp.initial_time] = unk[ctds.initial_time_index] : γ = unk[ctds.initial_time_index] 
  end
  if CTBase.__is_final_time_free(ctds.ocp)
    ctds.variable_dimension > 1 ? γ[ctds.ocp.final_time] = unk[ctds.final_time_index] : γ = unk[ctds.final_time_index]
  end

  for i ∈ 1:(ctds.grid_size_fine+1)
    if i ∈ 1:(ctds.grid_size_fine÷ctds.grid_size_coarse):(ctds.grid_size_fine+1)
      xᵢ = unk[rg(ctds.time_size + (k-1)*ctds.ocp.state_dimension + 1, ctds.time_size + k*ctds.ocp.state_dimension)]
      if i == 1
        c[ctds.time_size-1+1:ctds.time_size-1+ctds.ocp.state_dimension] = zeros(eltype(unk),ctds.ocp.state_dimension)
      else
        c[rg(ctds.time_size-1 + (k-1)*ctds.ocp.state_dimension + 1, ctds.time_size-1 + k*ctds.ocp.state_dimension)] = xᵢ-xᵣ
      end
      k += 1
    else
      xᵢ = xᵣ
    end

    tᵢ = times[i]

    if i < ctds.grid_size_fine+1
      tᵢ₊₁ = times[i+1]
      if ctds.fixed_time_step
        c[1] = unk[ctds.final_time_index]-unk[ctds.initial_time_index]
      else
        c[i] = tᵢ₊₁ - tᵢ
      end
      tspan = (tᵢ,tᵢ₊₁)
      uᵢ = ctds.u_fun(β,i)
      arg = [uᵢ;γ]
      if ctds.ocp.state_dimension > 1
        p_dynamics = ODEProblem(ctds.dyna!, xᵢ, tspan, arg)
      else
        p_dynamics = ODEProblem(ctds.dyna, xᵢ, tspan, arg)
      end
      sol_dynamics = OrdinaryDiffEq.solve(p_dynamics, Tsit5(), abstol=ctds.tol)
      xᵣ = sol_dynamics(tᵢ₊₁)


      if !isnothing(ctds.ocp.lagrange)
        p_lagrange = ODEProblem(ctds.f⁰, 0.0, tspan, [xᵢ;uᵢ;γ])
        sol_lagrange = OrdinaryDiffEq.solve(p_lagrange, Tsit5(), abstol=ctds.tol)
        lagrange_cost += sol_lagrange(tᵢ₊₁)
      end


    end

    offset = ctds.time_size-1 + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension + ctds.bc_dimension
    if xind_dim > 0
      c[offset + (i-1)*xind_dim + 1 : offset + i*xind_dim] = [xᵢ[ind] for ind ∈ xind]
    end
    
    offset = offset + (ctds.grid_size_fine+1)*xind_dim
    if uind_dim > 0
      c[offset + (i-1)*uind_dim + 1 : offset + i*uind_dim] = [uᵢ[ind] for ind ∈ uind]
    end
    
    offset = offset + (ctds.grid_size_fine+1)*uind_dim
    if ctds.nlc_dimension > 0
      c[offset + (i-1)*ctds.nlc_dimension + 1 : offset + i*ctds.nlc_dimension] = ctds.nlc(tᵢ,xᵢ,uᵢ,γ)
    end

  end

  offset = ctds.time_size-1 + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension
  if ctds.bc_dimension > 0
    c[offset+1:offset+ctds.bc_dimension] = ϕ(unk[rg(ctds.time_size + 1, ctds.time_size + ctds.ocp.state_dimension)], unk[rg(ctds.time_size + (ctds.grid_size_coarse)*ctds.ocp.state_dimension + 1, ctds.time_size + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension)], γ)
  end

  offset = offset + ctds.bc_dimension + (xind_dim + uind_dim + ctds.nlc_dimension)*(ctds.grid_size_fine+1)  
  if ctds.vc_dimension > 0
    c[offset+1:offset+ctds.vc_dimension] = θ(γ)
  end

  mayer_cost = 0
  if !isnothing(ctds.ocp.mayer)
    mayer_cost = ctds.ocp.mayer(unk[rg(ctds.time_size + 1, ctds.time_size + ctds.ocp.state_dimension)], unk[rg(ctds.time_size + (ctds.grid_size_coarse)*ctds.ocp.state_dimension + 1, ctds.time_size + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension)], γ)
  end

  c[end] = mayer_cost+lagrange_cost - unk[end]

  return c

end
