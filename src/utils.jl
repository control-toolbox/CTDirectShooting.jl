# --------------------------------------------------------------------------------------------------
# Utils for the transcription from ocp to descent problem

function parse_ocp_direct_shooting(ocp::OptimalControlModel)

    # parsing ocp
    dy = ocp.dynamics
    co = ocp.lagrange
    n = ocp.state_dimension
    m = ocp.control_dimension

    # initial_condition and final_constraint
    # x0 = ocp.initial_condition
    # cf = ocp.final_constraint
    x0 = nothing
    cf = nothing
    constraints = ocp.constraints
    for (_, c) ∈ constraints
        @match c begin
        (:initial, _, f, lb, ub) => begin
            if lb ≠ ub 
                error("direct shooting is implemented for problems with initial condition")
            else
                x0 = lb
            end
            end
        (:final, _, f, lb, ub) => begin
            if lb ≠ ub 
                error("direct shooting is implemented for problems with final equality constraint")
            else
                cf = x -> f(x) - lb
            end
            end
	    end # match
    end # for

    return dy, co, cf, x0, n, m

end

# forward integration of the state
"""
	model(x0, T, U, f)

TBW
"""
function model_primal_forward(x0, T, U, f)
    xₙ = x0
    X = [xₙ]
    for n in range(1, length(T) - 1)
        xₙ = f(T[n], xₙ, T[n+1], U[n])
        X = vcat(X, [xₙ]) # vcat gives a vector of vector
    end
    return xₙ, X
end

# backward integration of state and costate
"""
	adjoint(xₙ, pₙ, T, U, f)

TBW
"""
function model_adjoint_backward(xₙ, pₙ, T, U, f)
    X = [xₙ]
    P = [pₙ]
    for n in range(length(T), 2, step=-1)
        xₙ, pₙ = f(T[n], xₙ, pₙ, T[n-1], U[n-1])
        X = vcat([xₙ], X)
        P = vcat([pₙ], P)
    end
    return xₙ, pₙ, X, P
end