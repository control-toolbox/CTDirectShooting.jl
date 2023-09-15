using CTDirectShooting
using CTBase


t0 = 0
x0 = 0

@def ocp begin
    tf ∈ R, variable
    t ∈ [t0, tf], time
    x ∈ R, state
    u ∈ R, control
    x(t0) == x0, initial_con
    (x(tf) - tf) - 10 == 0, boundary_constraint
    ẋ(t) == u(t)
    ∫(0.5 * u(t) ^ 2) → min
end

sol = solve(ocp,grid_size_fine=30,grid_size_coarse=10,fixed_time_step=true)

println(sol)