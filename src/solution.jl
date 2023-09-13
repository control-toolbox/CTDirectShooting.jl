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

    function Solution(unk_sol::Vector{},prob::SimpleProblem,N::Int64,M::Int64)
        times = get_times(unk_sol,N)
        state = get_coarse_states(unk_sol,prob.state_dim,N,M)
        control = get_control(unk_sol,prob,N,M)
        variable = [0]
        objective = unk_sol[end]
        return new(times,state,control,variable,objective)
    end


    function Solution(unk_sol::Vector{}, ocp::OptimalControlModel, parameter_dimension, N::Int64, M::Int64)
        times = get_times(unk_sol,N)
        state = get_coarse_states(unk_sol,ocp.state_dimension,N,M)
        control = get_control(unk_sol,ocp,N,M)
        variable = get_variable(unk_sol,ocp,parameter_dimension,N,M)
        objective = unk_sol[end]
        return new(times,state,control,variable,objective)
    end


end