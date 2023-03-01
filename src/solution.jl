
# ------------------------------------------------------------------------------------
# Direct shooting solution
#
struct DirectShootingSolution <: AbstractOptimalControlSolution
    T::TimesDisc # the times
    X::States # the states at the times T
    U::Controls # the controls at T
    P::Adjoints # the adjoint at T
    objective::MyNumber
    state_dimension::Dimension # the dimension of the state
    control_dimension::Dimension # the dimension of the control
    stopping::Symbol # the stopping criterion
    message::String # the message corresponding to the stopping criterion
    success::Bool # whether or not the method has finished successfully: CN1, stagnation vs iterations max
    iterations::Integer # the number of iterations
end

function DirectShootingSolution(sol::CTOptimization.UnconstrainedSolution,
    ocp::OptimalControlModel, grid::TimesDisc, penalty_constraint::Real)

    # 
    VFN = VectorField{:nonautonomous}

    # parsing ocp
    dy, co, cf, x0, n, m = parse_ocp_direct_shooting(ocp)

    # control solution
    U⁺ = vec2vec(sol.x, m)

    # Jacobian of the constraints
    Jcf(x) = ctjacobian(cf, x)

    # penalty term for final constraints
    αₚ = penalty_constraint

    # state flow
    f = Flow(VFN(dy))

    # augmented state flow
    fa = Flow(VFN((t, x, u) -> [dy(t, x[1:end-1], u)[:]; co(t, x[1:end-1], u)]))

    # flow for state-adjoint
    p⁰ = -1.0
    H(t, x, p, u) = p⁰ * co(t, x, u) + p' * dy(t, x, u)
    fh = Flow(Hamiltonian{:nonautonomous}(H))

    # get state and adjoint
    T = grid
    xₙ, _ = model_primal_forward(x0, T, U⁺, f)
    pₙ = p⁰ * αₚ * transpose(Jcf(xₙ)) * cf(xₙ)
    _, _, X⁺, P⁺ = model_adjoint_backward(xₙ, pₙ, T, U⁺, fh)

    # function J, that we minimize
    function J(U::Controls)
        # via augmented system
        xₙ, X = model_primal_forward([x0[:]; 0.0], T, U, fa)
        cost = xₙ[end] + 0.5 * αₚ * norm(cf(xₙ[1:end-1]))^2
        return cost
    end
    objective = J(U⁺)

    dssol = DirectShootingSolution(T, X⁺, U⁺, P⁺, objective, n, m, 
        sol.stopping, sol.message, sol.success, sol.iterations)
   
    return _OptimalControlSolution(ocp, dssol)

end

function _OptimalControlSolution(ocp::OptimalControlModel, dssol::DirectShootingSolution)
    x = ctinterpolate(dssol.T, dssol.X) # je ne peux pas donner directement la sortie de ctinterpolate car ce n'est pas une Function
    u = ctinterpolate(dssol.T[1:end-1], dssol.U)
    p = ctinterpolate(dssol.T, dssol.P)
    sol = OptimalControlSolution()
    sol.state_dimension = dssol.state_dimension
    sol.control_dimension = dssol.control_dimension
    sol.times = dssol.T
    sol.state = t -> x(t)
    sol.state_labels = ocp.state_labels # update CTBase to have a getter
    sol.adjoint = t -> p(t)
    sol.control = t -> u(t)
    sol.control_labels = ocp.control_labels
    sol.objective = dssol.objective
    sol.iterations = dssol.iterations
    sol.stopping = dssol.stopping
    sol.message = dssol.message
    sol.success = dssol.success
    return sol
end