"""
$(TYPEDEF)

**Fields**

$(TYPEDFIELDS)

"""
struct Solution
    times::Vector{}
    state::Vector{}
    control::Vector{}
    variable::Vector{}
    objective::Real

    function Solution(unk_sol::Vector{},prob::SimpleProblem,N::Int64,M::Int64)
        times = get_times(unk_sol,N)
        state = get_coarse_states(unk_sol,prob.state_dim,N,M)
        control = get_control(unk_sol,prob,N,M)
        variable = [0]
        objective = unk_sol[end]
        return new(times,state,control,variable,objective)
    end

end