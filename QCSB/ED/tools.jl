function pad(M::AbstractMatrix, before::Int, after::Int)
    return kron(sparse(I, 2^before, 2^before), M, sparse(I, 2^after, 2^after))
end