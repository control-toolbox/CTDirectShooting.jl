"""
$(TYPEDSIGNATURES)

OptimalControlSolution generation
"""
function _OptimalControlSolution(unk_sol, ctds)
    times = get_times(unk_sol,ctds)
    controls = reshape_vector(get_control(unk_sol,ctds),ctds.ocp.control_dimension)
    states = reshape_vector(get_states(unk_sol,ctds),ctds.ocp.state_dimension)
    variable = is_variable_dependent(ctds.ocp) ? get_variable(unk_sol,ctds) : 0.0
    x = ctinterpolate(times, states)
    u = ctinterpolate(times, controls)
    p = ctinterpolate(times, states)
    v = variable
    sol = OptimalControlSolution() # +++ constructor with ocp as argument ?
    copy!(sol, ctds.ocp)
    sol.times      = times
    sol.state      = (sol.state_dimension==1)    ? deepcopy(t -> x(t)[1]) : deepcopy(t -> x(t)) # scalar output if dim=1
    sol.costate    = (sol.state_dimension==1)    ? deepcopy(t -> p(t)[1]) : deepcopy(t -> p(t)) # scalar output if dim=1
    sol.control    = (sol.control_dimension==1)  ? deepcopy(t -> u(t)[1]) : deepcopy(t -> u(t)) # scalar output if dim=1
    sol.variable   = (sol.variable_dimension==1) ? v[1] : v # scalar output if dim=1
    sol.objective  = unk_sol[end]
    sol.iterations = 0 #ctd.NLP_iterations todo
    sol.stopping   = :dummy
    sol.message    = "no message"
    sol.success    = false #

    return sol
end