function reduced(state::DiagonalStateED, A::AbstractVector{Int})
    ρ = state.vec
    L = Int(log2(length(ρ)))

    A = [L-site+1 for site in A]  # Reverse order for ITensors.jl

    # Convert state vector to multidimensional array
    state_tensor = reshape(state, fill(2, L)...)
    
    # Get complement of sites
    B = setdiff(1:L, A)
    
    # Sum over all other sites
    result = dropdims(sum(state_tensor, dims=B), dims=Tuple(B))

    # Convert back to vector
    return DiagonalStateED(vec(result))
end

function partial_trace(state::DiagonalStateED, B::AbstractVector)
    A = setdiff(1:Int(log2(length(state.vec))), B)
    return reduced(state, A)
end

function entropy(state::DiagonalStateED)
    ps = state.vec
    return -sum([p*log2(p) for p in ps if real(p) > 0])
end

# well suited for generalization to be implementation agnostic
# just need general length() function 

function mutual_info(state::DiagonalStateED, A::AbstractVector, B::AbstractVector)
    C = [i for i in 1:Int(log2(length(state.vec))) if !(i in A) && !(i in B)]
    SA = entropy(partial_trace(state, union(B,C)))
    SB = entropy(partial_trace(state, union(A,C)))
    SAB = entropy(partial_trace(state, C))
    return SA + SB - SAB
end