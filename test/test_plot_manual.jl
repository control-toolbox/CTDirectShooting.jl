using CTDirectShooting
using CTProblemLibrary
using CTBase

prob = Problem(:integrator, :dim2, :energy); ocp = prob.model
sol = solve(ocp) # print problem
plot(sol)