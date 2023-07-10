# --------------------------------------------------------------------------------------------------
# make a simple ocp model
# direct simple shooting
struct SimpleProblem
    g::Function # mayer cost function
    f⁰::Function # langrange cost function
    f::Function # dynamic function
    dmin::Vector{} # lower bound box constraints (state, control)
    dmax::Vector{} # upper bound box constraints (state, control)
    cmin::Vector{} # lower bound constraints (state, control)
    c::Function # constraints function (state, control, var)
    cmax::Vector{} # upper bound constraints (state, control)
    bmin::Vector{} # lower bound boundary conditions (t0, x(t0), tf, x(tf))
    b::Function # boundary conditions function (t0, x(t0), tf, x(tf), var)
    bmax::Vector{} # upper bound boundary conditions (t0, x(t0), tf, x(tf))

    function SimpleProblem(g::Function, f⁰::Function, f::Function, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Function, cmax::Vector{}, bmin::Vector{}, b::Function, bmax::Vector{})

        
        @assert size(dmin) == size(dmax) 
        #@assert typeof(dmin) == typeof(dmax)   
        @assert size(cmin) == size(cmax) 
        #@assert typeof(cmin) == typeof(cmax) 
        @assert size(bmin) == size(bmax) 
        #@assert typeof(bmin) == typeof(bmax) 
    
        return new(g::Function, f⁰::Function, f::Function, dmin::Vector{}, dmax::Vector{}, cmin::Vector{}, c::Function, cmax::Vector{}, bmin::Vector{}, b::Function, bmax::Vector{})
    
    end
    
end
