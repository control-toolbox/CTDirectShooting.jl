function solve(prob::SimpleProblem, N::Int64=9, M::Int64=N÷3)
    
    # checking
    _check_grids(N,M)

    boundary_freedom = BoundaryFreedom(prob)

    u(β,i) = prob.l(β[i]) #picewise constant  
    
    obj_fun(unk) = unk[N+1] #
    dim_unk = size(I,1) + 2*size(J,1) + size(I,1)-1 # N+1 + 2*M + 1*N
    unk0 = fill(0.5,dim_unk)
    l_var = -Inf*ones(dim_unk)
    u_var = Inf*ones(dim_unk)
    for i ∈ 1:dim_unk # variable bounds
        if i ≤ N+1
            l_var[i] = 0
            u_var[i] = Inf
        end
    end

    return sol
end