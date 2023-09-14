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


"""
$(TYPEDSIGNATURES)

return controls

"""
function get_variable(unk,ocp::OptimalControlModel,parameter_dimension,N,M,time_size) 
    variable = unk[rg(time_size + M*ocp.state_dimension + N*parameter_dimension + 1,time_size + M*ocp.state_dimension + N*parameter_dimension + ocp.variable_dimension)]
    return variable
end


# return a view of unknowns state on coarse grid
"""
$(TYPEDSIGNATURES)

return a view of unknowns state on coarse grid 

"""
get_coarse_states(unk,dim,N,M,time_size) = view(unk,time_size+1:time_size+M*dim)

# return a view of unknowns parameters β 
"""
$(TYPEDSIGNATURES)

return a view of unknowns parameters β 

"""
get_parameters(unk,dim1,dim2,N,M, time_size) = view(unk,time_size+dim1*M+1:time_size+dim1*M+dim2*N)

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

function OCP_To_SimpleProblem(ocp::OptimalControlModel)
    zero_fun(;kwargs) = zero(eltype(kwargs))
    g = isnothing(ocp.mayer) ? zero_fun : ocp.mayer
    f⁰ = isnothing(ocp.lagrange) ?  zero_fun : ocp.langrange
    f = ocp.dynamics
    function f!(dx, x, u, t)
      for i in range(1,size(x,1))
        dx[i] = f(t, x, u)[i]
      end
    end
    state_dim = ocp.state_dimension
    control_dim = ocp.control_dimension
    (ξl, ξ, ξu), (ηl, η, ηu), (ψl, ψ, ψu), (ϕl, ϕ, ϕu), (θl, θ, θu), (ul, uind, uu), (xl, xind, xu), (vl, vind, vu) = nlp_constraints(ocp)

    b(t0,tf,x0,xf) = [t0,tf,ϕ(x0,xf)]
    bmin = fill(-Inf,6)
    bmax = fill(Inf,6)
    CTBase.__is_initial_time_free(ocp) ? nothing : bmin[1] = ocp.initial_time
    CTBase.__is_final_time_free(ocp) ? nothing : bmin[2] = ocp.final_time
    bmin[3:end] = ϕl
    CTBase.__is_initial_time_free(ocp) ? nothing : bmax[1] = ocp.initial_time
    CTBase.__is_final_time_free(ocp) ? nothing : bmax[2] = ocp.final_time
    bmax[3:end] = ϕu

    d = [xind;uind]
    dmin = [xl;ul]
    dmax = [xu;uu]

    c = [η,ξ,ψ]
    cmin = [ηl,ξl,ψl]
    cmax = [ηu,ξu,ψu]

    l = β -> β
    # dmin = Vector{eltype(ocp.initial_time)}(-Inf,4)
    # dmin[]

    return SimpleProblem(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any)
  end