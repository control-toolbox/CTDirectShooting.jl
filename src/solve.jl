# availble methods by order of preference: from top to bottom
algorithmes = ()
algorithmes = add(algorithmes, (:adnlp, :ipopt))

"""
$(TYPEDSIGNATURES)

Return the list of available methods to solve the optimal control problem.
"""
function available_methods()::Tuple{Tuple{Vararg{Symbol}}}
    return algorithmes
end


"""
$(TYPEDSIGNATURES)

Solve a SimpleProblem with discretisation and return the solution : 
```(tᵢ)[0:N=I], (xⱼ)[J], (βᵢ)[0:N-1], var, obj```
"""
function solve(ocp::OptimalControlModel,
  description...;
  display::Bool=__display(),
  grid_size_fine::Integer=__grid_size_fine,
  grid_size_coarse::Integer=__grid_size_coarse(grid_size_fine),
  tol::Float64=1e-9,
  fixed_time_step::Bool=true,
  print_level::Integer=__print_level_ipopt(),
  mu_strategy::String=__mu_strategy_ipopt(),
  kwargs...)
  
  # check grids dimensions
  _check_grids(grid_size_fine,grid_size_coarse)

  # get full description from partial
  # throw error if description is not valid
  # thus, no else is needed below
  method = getFullDescription(description, available_methods())
  
  if :adnlp in method
    ctds = CTDirectShooting_data(ocp,grid_size_fine,grid_size_coarse,fixed_time_step,init,tol)
    unk0 = _initial_guess(ctds,0.5)
    #unk0 = _init_vector(ctds,0.5,ctds.time_size)
    l_var, u_var = variable_bounds(ctds)
    lb, ub = constraint_bounds(ctds)
    
    nlp = ADNLPModel!(unk -> obj_fun(unk,ctds), unk0, l_var, u_var, (c,unk) -> con_fun(c,unk,ctds), lb, ub)
  end

  # solve
  if :ipopt in method
    # https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl/blob/main/src/NLPModelsIpopt.jl#L119
    # callback: https://github.com/jump-dev/Ipopt.jl#solver-specific-callback
    # sb="yes": remove ipopt header +++ make that default
    print_level = display ?  print_level : 0
    solver_ad = IpoptSolver(nlp)
    ipopt_solution_ad = NLPModelsIpopt.solve!(solver_ad, nlp, tol=tol, mu_strategy=mu_strategy, sb="yes", print_level=print_level)
end
  
  sol = _OptimalControlSolution(ipopt_solution_ad.solution, ctds)

  return sol

end