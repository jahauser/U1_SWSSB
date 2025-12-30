# "doubled" MPS corresponding to the identity operator on a single site
# in practice, this is a bell pair (without normalization)
function id_mps(s1::Index, s2::Index)
    return MPS([1; 0; 0; 1], [s1, s2])
end

# "doubled" MPS corresponding to the identity operator on N sites
# in practice, this looks like the tensor product of many bell pairs
function bell(sites::Vector{Index{Int}})
    N = length(sites)
    tensors = ITensor[]
    for i in 1:2:N
        a, b = id_mps(sites[i], sites[i+1])
        push!(tensors, a, b)
    end 
    return MPS(tensors)
end

# "doubled" MPS corresponding to operator M at site pos and the identity everywhere else
function M_bra(sites::Vector{Index{Int}}, M::AbstractMatrix, pos::Int; refs=0)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - refs

    bra = bell(sites)
    bra = apply(op(M, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...), bra)
    return bra
end

# the trace of the density matrix encoded in a "doubled" MPS
function doubledtrace(ρ::MPS)
    return inner(bell(siteinds(ρ)), ρ)
end