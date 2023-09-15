using CTDirectShooting
using CTBase


t0=0
tf=1
x0=[-1, 0]
xf=[0, 0]
A = [ 0 1
    0 0 ]
B = [ 0
1 ]

@def ocp begin
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control
    x(t0) == [-1, 0], initial_con
    x(tf) == [0, 0], final_con
    ẋ(t) == A * x(t) + B * u(t)
    ∫(0.5 * u(t) ^ 2) → min
end

sol = solve(ocp,grid_size_fine=30,grid_size_coarse=3,fixed_time_step=true)

println(sol)