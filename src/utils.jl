# range same as : but return i if i:i
"""
$(TYPEDSIGNATURES)

range same as : but return i if i:i

"""
rg(i,j) = i == j ? i : i:j

# return the element i of size l in a vector
"""
$(TYPEDSIGNATURES)

return the element i of size l in a vector

"""
get_element(v,i,l) = view(v,(l*(i-1) + 1):l*i)

# return a vector of uniform times between t0 and tf 
"""
$(TYPEDSIGNATURES)

return a vector of uniform times between t0 and tf 

"""
function get_times_uniform(N,t0,tf) # case t0 fixed
    return [(t0*(N+1-i)+tf*(i-1))/N for i ∈ 1:N+1]
end

# return a view of unknowns times 
"""
$(TYPEDSIGNATURES)

return a view of unknowns times 

"""
get_times(unk,N,fixed_time_step) = fixed_time_step ? get_times_uniform(N,unk[1],unk[2]) : view(unk,1:N+1)

# return a view of unknowns times 
"""
$(TYPEDSIGNATURES)

return a view of unknowns times 

"""
get_times(unk,ctds) = ctds.fixed_time_step ? get_times_uniform(ctds.grid_size_fine,unk[1],unk[2]) : view(unk,1:ctds.grid_size_fine+1)

# return controls 
# """
# $(TYPEDSIGNATURES)

# return controls

# """
# function get_control(unk,prob::SimpleProblem,N,M) 
#     control = [prob.l(view(unk,N+1+M*prob.state_dim+i:N+1+M*prob.state_dim+i+prob.control_dim-1)) for i in 1:prob.control_dim:N]
#     return control
# end

"""
$(TYPEDSIGNATURES)

return controls

"""
function get_control(unk,ocp::OptimalControlModel,N,M,time_size) 
    control = [view(unk,time_size+M*ocp.state_dimension+i:time_size+M*ocp.state_dimension+i+ocp.control_dimension-1) for i in 1:ocp.control_dimension:N]
    return control
end
function get_control(unk,ctds) # only constant control for now
    β = get_parameters(unk,ctds)
    return β
end


"""
$(TYPEDSIGNATURES)

return variables

"""
function get_variable(unk,ocp::OptimalControlModel,parameter_dimension,N,M,time_size) 
    variable = unk[rg(time_size + M*ocp.state_dimension + N*parameter_dimension + 1,time_size + M*ocp.state_dimension + N*parameter_dimension + ocp.variable_dimension)]
    return variable
end
function get_variable(unk,ctds) 
    variable = unk[rg(ctds.time_size + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension + ctds.grid_size_fine*ctds.parameter_dimension + 1,ctds.time_size + (ctds.grid_size_coarse+1)*ctds.ocp.state_dimension + ctds.grid_size_fine*ctds.parameter_dimension + ctds.ocp.variable_dimension)]
    return variable
end


# return a view of unknowns state on coarse grid
"""
$(TYPEDSIGNATURES)

return a view of unknowns state on coarse grid 

"""
get_coarse_states(unk,dim,N,M,time_size) = view(unk,time_size+1:time_size+M*dim)
get_coarse_states(unk,ctds) = view(unk,ctds.time_size+1:ctds.time_size+(ctds.grid_size_coarse+1)*ctds.ocp.state_dimension)

"""
$(TYPEDSIGNATURES)

return a vector of states on fine grid 

"""
function get_states(unk,ctds)
    states = zeros(0)
    times = get_times(unk,ctds)
    
    k=1
    xᵣ = ctds.ocp.state_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ctds.ocp.state_dimension)
    uᵢ = ctds.ocp.control_dimension==1 ? zero(eltype(unk)) : zeros(eltype(unk),ctds.ocp.control_dimension)
    β = get_parameters(unk,ctds)
    γ = ctds.variable_dimension > 0 ? get_variable(unk,ctds) : Vector{eltype(unk)}()

    for i ∈ 1:(ctds.grid_size_fine+1)
        if i ∈ 1:(ctds.grid_size_fine÷ctds.grid_size_coarse):(ctds.grid_size_fine+1)
          xᵢ = unk[rg(ctds.time_size + (k-1)*ctds.ocp.state_dimension + 1, ctds.time_size + k*ctds.ocp.state_dimension)]
          k += 1
        else
          xᵢ = xᵣ
        end
    
        tᵢ = times[i]
    
        if i < ctds.grid_size_fine+1
          tᵢ₊₁ = times[i+1]
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
        end
        append!(states,xᵢ)
    end    
    return states
end

# return a view of unknowns parameters β 
"""
$(TYPEDSIGNATURES)

return a view of unknowns parameters β 

"""
get_parameters(unk,dim1,dim2,N,M, time_size) = view(unk,time_size+dim1*M+1:time_size+dim1*M+dim2*N)
"""
$(TYPEDSIGNATURES)

return a view of unknowns parameters β 

"""
get_parameters(unk,ctds) = view(unk,ctds.time_size+ctds.ocp.state_dimension*(ctds.grid_size_coarse+1)+1:ctds.time_size+ctds.ocp.state_dimension*(ctds.grid_size_coarse+1)+ctds.parameter_dimension*ctds.grid_size_fine)

# return the index associated to a time
"""
$(TYPEDSIGNATURES)

return the index associated to a time

"""
function time_to_index(times,t)
    for i in size(times,1)-1
        if times[i] ≤ t < times[i+1]
            return i
        else t ≤ times[i+1]
            return i+1
        end
    end
end

"""
$(TYPEDSIGNATURES)

return the initial guess

"""
function _initial_guess(ctds,val::T) where T<:Real
    
    vec = fill(val, ctds.unk_dim) 
    if is_variable_dependent(ctds.ocp)
        γ = get_variable(vec,ctds)
    end
    if CTBase.__is_initial_time_free(ctds.ocp)
        t0 = γ[ctds.ocp.initial_time] 
    else
        t0 = ctds.ocp.initial_time
    end
    if CTBase.__is_final_time_free(ctds.ocp)
        tf = γ[ctds.ocp.final_time]
    else
        tf = ctds.ocp.final_time
    end
    for i ∈ 1:(ctds.time_size)
        vec[i] = (t0*(ctds.time_size-i)+tf*(i-1))/(ctds.time_size-1)
        # add other inits
    end
    return vec 
end