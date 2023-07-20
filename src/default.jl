# --------------------------------------------------------------------------------------------------
# Direct shooting method - default values

"""
$(TYPEDSIGNATURES)

Used to set the default value of the grid size for the direct shooting method.
The default value is `201`.
"""
__grid_size_direct_shooting() = 201

"""
$(TYPEDSIGNATURES)

Used to set the default value of the penalty term in front of the final constraint for the direct shooting method.
The default value is `1e4`.
"""
__penalty_term_direct_shooting() = 1e4

"""
$(TYPEDSIGNATURES)

Used to set the default value of the maximal number of iterations for the direct shooting method.
The default value is `100`.
"""
__max_iter_direct_shooting() = 100

"""
$(TYPEDSIGNATURES)

Used to set the default value of the absolute tolerance for the stopping criterion for the direct shooting method.
The default value is `10 * eps()`.
"""
__abs_tol_direct_shooting() = 10 * eps()

"""
$(TYPEDSIGNATURES)

Used to set the default value of the optimality relative tolerance for the stopping criterion for the direct shooting method.
The default value is `1e-8`.
"""
__opt_tol_direct_shooting() = 1e-8

"""
$(TYPEDSIGNATURES)

Used to set the default value of the step stagnation relative tolerance for the stopping criterion for the direct shooting method.
The default value is `1e-8`.
"""
__stagnation_tol_direct_shooting() = 1e-8
