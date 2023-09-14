function _init_vector(N::Int64, val::T, dim::Int64, t0::Int64=0, tf::Int64=t0+1) where T <: Real 
    vec = fill(val, dim) 
    for i ∈ 1:N+1
        vec[i] = (t0*(N+1-i)+tf*(i-1))/N
    end
    return vec 
end