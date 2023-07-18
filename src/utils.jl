# range same as : but return i if i:i
rg(i,j) = i == j ? i : i:j

# return a view of unknowns times 
get_times(unk) = view(unk,1:N+1)

# return a vector of uniform times between t0 and tf 
function get_times_uniform(prob) # case t0 fixed
    if boundary_freedom.tf 
        t0 = prob.bmin[1]
        tf = t0 + 1
    else
        t0, tf = prob.bmin[1], prob.bmin[3] 
    end
    return [t0*(N-i)/N + i/N*tf for i ∈ 0:N]
end