using CTDirectShooting

# the problem (double integrator time control constraint)
t0 = 0
x0 = [-1, 0]
xf = [0, 0]
γ = 1

g(t0, x0, tf, xf) = tf
f⁰(t, x, u) = 0
f(t, x, u) = [x(t)[2],u(t)]
dmin = [-Inf, -Inf, -γ]
dmax = [Inf, Inf, γ]
cmin = []
c(t, x, u) = nothing
cmax = []
bmin = [0, -1, 0, 0, 0, 0]
b(t0, x0, tf, xf) = [t0, x0, tf, xf]
bmax = [0, -1, 0, Inf, 0, 0]

prob = SimpleProblem(g, f⁰, f, dmin, dmax, cmin, c, cmax, bmin, b, bmax)

# solving 

tf = t0+1
N = 9
I = [i for i ∈ 0:N]
T = range(t0, tf, N)
β = zeros(N) .+ 0.1
function u(t)
    index = 1
    while (t >= T[index] && index < N)
        index += 1
    end
    return(β[index])
end
U = [u(t) for t in T]
J = I[1:3:size(I,1)]
XJ = [(((lastindex(J)-j)/size(J,1))*x0 + j/size(J,1)*xf) for j in eachindex(J)] 
F = tf












nothing