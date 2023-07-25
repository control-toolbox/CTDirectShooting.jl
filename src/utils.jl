# range same as : but return i if i:i
rg(i,j) = i == j ? i : i:j

# return the element i of size l in a vector
get_element(v,i,l) = view(v,(l*(i-1) + 1):l*i)

# return a vector of uniform times between t0 and tf 
function get_times_uniform(prob,N) # case t0 fixed
    boundary_freedom = BoundaryFreedom(prob)
    if boundary_freedom.tf 
        t0 = prob.bmin[1]
        tf = t0 + 1
    else
        t0, tf = prob.bmin[1], prob.bmin[3] 
    end
    return [t0*(N-i)/N + i/N*tf for i ∈ 0:N]
end

# return a view of unknowns times 
get_times(unk,N) = view(unk,1:N+1)

# return a view of unknowns state on coarse grid
get_coarse_states(unk,dim,N,M) = view(unk,N+1+1:N+1+M*dim)

# return a view of unknowns parameters β 
get_parameters(unk,dim1,dim2,N,M) = view(unk,N+1+dim1*M+1:N+1+dim1*M+dim2*N)