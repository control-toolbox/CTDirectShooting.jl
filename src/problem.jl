# --------------------------------------------------------------------------------------------------
# make a simple ocp model
# direct simple shooting
"""
$(TYPEDEF)

**Fields**

$(TYPEDFIELDS)

"""
struct SimpleProblem
  g::Any # mayer cost function
  f⁰::Any # langrange cost function
  f::Any # dynamic function
  f!::Any # in-place dynamic function
  dmin::Vector{} # lower bound box constraints (state, control)
  dmax::Vector{} # upper bound box constraints (state, control)
  cmin::Vector{} # lower bound constraints (state, control)
  c::Any # constraints function (state, control, var)
  cmax::Vector{} # upper bound constraints (state, control)
  bmin::Vector{} # lower bound boundary conditions (t0, tf, x(t0), x(tf))
  b::Any # boundary conditions function (t0, tf, x(t0), x(tf), var)
  bmax::Vector{} # upper bound boundary conditions (t0, tf, x(t0), x(tf))
  state_dim::Int64
  control_dim::Int64
  l::Any # function to control u wrt β

  function SimpleProblem(g::Any, f⁰::Any, f::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any=(β -> β))


    @assert size(dmin) == size(dmax)
    #@assert typeof(dmin) == typeof(dmax)   
    @assert size(cmin) == size(cmax)
    #@assert typeof(cmin) == typeof(cmax) 
    @assert size(bmin) == size(bmax)
    #@assert typeof(bmin) == typeof(bmax) 

    #create the inplace dynamic function
    function f!(dx, x, u, t)
      # for i in range(1,size(x,1))
      #   dx[i] = f(t, x, u)[i]
      # end
      dx[1] = x[2]
      dx[2] = u
    end

    return new(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any)
  end

  function SimpleProblem(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any=(β -> β))


    @assert size(dmin) == size(dmax)
    #@assert typeof(dmin) == typeof(dmax)   
    @assert size(cmin) == size(cmax)
    #@assert typeof(cmin) == typeof(cmax) 
    @assert size(bmin) == size(bmax)
    #@assert typeof(bmin) == typeof(bmax) 

    return new(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any)
  end

  function SimpleProblem(ocp::OptimalControlModel)
    g = ocp.mayer
    f⁰ = ocp.langrange
    f = ocp.dynamic
    function f!(dx, x, u, t)
      for i in range(1,size(x,1))
        dx[i] = f(t, x, u)[i]
      end
    end

    (ξl, ξ, ξu), (ηl, η, ηu), (ψl, ψ, ψu), (ϕl, ϕ, ϕu), (θl, θ, θu), (ul, uind, uu), (xl, xind, xu), (vl, vind, vu) = nlp_constraints(ocp)

    b(t0,tf,x0,xf) = [t0,tf,ϕ(x0,xf)]
    bmin = fill(-Inf,4)
    bmax = fill(Inf,4)
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

    # dmin = Vector{eltype(ocp.initial_time)}(-Inf,4)
    # dmin[]

    return new(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any)
  end
    
    

  #   return new(g::Any, f⁰::Any, f::Any, f!::Any, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Any, cmax::Vector{}, bmin::Vector{}, b::Any, bmax::Vector{}, state_dim::Int64, control_dim::Int64, l::Any)
  # end

end

"""
$(TYPEDEF)

**Fields**

$(TYPEDFIELDS)

"""
struct BoundaryFreedom
  t0::Bool
  tf::Bool
  x0::Bool
  xf::Bool
  function BoundaryFreedom(prob::SimpleProblem)
    bv = prob.bmin .!= prob.bmax
    return new(bv[1], bv[2], bv[3] && bv[4], bv[5] && bv[6])
  end
end
