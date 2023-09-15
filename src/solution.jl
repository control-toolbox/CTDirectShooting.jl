"""
$(TYPEDEF)

**Fields**

$(TYPEDFIELDS)

"""
struct Solution
    times::Vector{}
    state::Vector{}
    control::Vector{}
    variable
    objective::Real


    function Solution(unk_sol::Vector{}, ocp::OptimalControlModel, parameter_dimension, N::Int64, M::Int64, time_size)
        times = get_times(unk_sol,N,time_size==2)
        state = get_coarse_states(unk_sol,ocp.state_dimension,N,M,time_size)
        control = get_control(unk_sol,ocp,N,M,time_size)
        variable = is_variable_dependent(ocp) ? get_variable(unk_sol,ocp,parameter_dimension,N,M, time_size) : 0
        objective = unk_sol[end]
        return new(times,state,control,variable,objective)
    end


end

function _OptimalControlSolution(unk_sol, ctds)
    times = get_times(unk_sol,ctds)
    controls = get_control(unk_sol,ctds)
    states = get_states(unk_sol,ctds)
    x = ctinterpolate(times, states)
    u = ctinterpolate(times, controls)
    p = ctinterpolate(times, states)
    sol = OptimalControlSolution() # +++ constructor with ocp as argument ?
    copy!(sol, ocp)
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
end