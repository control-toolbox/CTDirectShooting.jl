function _init_vector(N::Int, val::T, dim::Int64) where T <: Real # change to accept custom t0/tf, here t0 = 0 and tf = 1
    vec = fill(val, dim) 
    for i ∈ 1:N+1
        vec[i] = 1.0*(i-1)/N
    end
    return Vector{Real}(vec) 
end