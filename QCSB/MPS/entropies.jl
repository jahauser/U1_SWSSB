function Svn(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    _, Λ, _ = svd(ψ[b], linkind(ψ, b))

    ps = [Λ[i,i]^2 for i in 1:dim(Λ,1) if Λ[i,i] > 0]
    return -sum([p*log2(p) for p in ps])
end

function reduced(state::DiagonalStateMPS, A::AbstractVector{Int})
    B = setdiff(1:length(siteinds(state.mps)), A)
    return partial_trace(state, B)
end

function partial_trace(state::DiagonalStateMPS, B::AbstractVector)
    ρ = state.mps
    sites = siteinds(ρ)

    ρA = Vector{ITensor}()
    B = sort(B)

    i = 1
    left = ITensor(1)
    while i <= length(B) && B[i] == i
        left = left * ITensor([1, 1], sites[i]) * ρ[i]
        i += 1
    end

    push!(ρA, ρ[i])
    ρA[1] *= left 

    for site in i+1:length(sites)
        if site in B
            ρA[end] *= ITensor([1, 1], sites[site]) * ρ[site]
        else
            push!(ρA, ρ[site])
        end
    end

    return DiagonalStateMPS(MPS(ρA))
end

function entropy(state::DiagonalStateMPS)
    ρ = dense(state)
    ps = diag(ρ)
    return -sum([p*log2(p) for p in ps if real(p) > 0])
end

function mutual_info(state::DiagonalStateMPS, A::AbstractVector, B::AbstractVector)
    C = [i for i in 1:length(siteinds(state.mps)) if !(i in A) && !(i in B)]
    SA = entropy(partial_trace(state, union(B,C)))
    SB = entropy(partial_trace(state, union(A,C)))
    SAB = entropy(partial_trace(state, C))
    return SA + SB - SAB
end