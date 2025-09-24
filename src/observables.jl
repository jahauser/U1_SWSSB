function entanglemententropy(ψ::MPS; bond = nothing)
    # number of qubits
    L = length(ψ)

    # make sure the state is normalized
    ψ = normalize!(ψ)

    # if bond = nothing, take the center bond
    bond = (isnothing(bond) ? L÷2 : bond)
    @assert (bond < L)

    # gauge the MPS
    orthogonalize!(ψ, bond)

    # get singular values
    row_inds = (bond > 1 ? (linkind(ψ, bond-1), siteind(ψ,bond)) : siteind(ψ,bond))
    _,s,_ = svd(ψ[bond],row_inds)

    # Compute Von Neumann Entropy S = -Tr(ρ log(ρ))
    S = 0.0
    for n in 1:dim(s, 1)
        λ = s[n,n]^2
        S -= λ * log(2, λ + 1e-20)
    end
    return S
end

function χ2_corr(Φ::MPS)
    sites = siteinds(Φ)
    L = length(Φ)
    Φ /= inner(MPS(sites, i -> "+"), Φ)*2^(L/2)
    return correlation_matrix(Φ, "a", "a†")
end

function χ2(Φ::MPS)
    L = length(Φ)

    corrs = χ2_corr(Φ)
    corrs[diagind(corrs)] .= 0
    χ = sum(corrs)/L + 0.25
    return χ
end

function entropy(state::DiagonalState, sites::AbstractVector, samples::Int; cutoff=1E-8, maxdim=200)
    entropy1 = 0.0
    entropy2 = 0.0

    for i in 1:samples
        state_copy = deepcopy(state) # I *think* this is redundant but not 100% sure, should figure out
        _, ms, vals = measure(state_copy, PauliZ, 1.0, sites; cutoff=cutoff, maxdim=maxdim)

        logprob = sum([!ms[i] ? log( (1+vals[i])/2 ) : log( (1-vals[i])/2 ) for i in eachindex(ms)])

        entropy1 += -logprob
        entropy2 += logprob^2
    end

    return entropy1/samples, entropy2/samples - (entropy1/samples)^2
end

function entropy(ρ::AbstractVector)
    return sum([real(val) > 0.0 ? -val * log(val) : 0.0 for val in ρ])
end

function reduced(state::DiagonalState, i::Int, j::Int)
    ψ = state.mps
    sites = siteinds(ψ)

    if i > 1
        left = ITensor([1, 1], sites[1]) * ψ[1]
    else
        left = ITensor(1)
    end
    for pos in 2:i-1
        left = left * ITensor([1, 1], sites[pos]) * ψ[pos]
    end

    if j < length(sites)
        right = ITensor([1, 1], sites[end]) * ψ[end]
    else
        right = ITensor(1)
    end
    for pos in length(sites)-1:-1:j+1
        right = right * ITensor([1, 1], sites[pos]) * ψ[pos]
    end

    return to_vector(left * contract(ψ[i:j]) * right)
end

function reduced(state::AbstractVector, sites::AbstractVector{Int})
    L = Int(log2(length(state)))

    sites = [L-site+1 for site in sites]  # Reverse order for ITensors.jl

    # Convert state vector to multidimensional array
    state_tensor = reshape(state, fill(2, L)...)
    
    # Get complement of sites
    other_sites = setdiff(1:L, sites)
    
    # Sum over all other sites
    result = dropdims(sum(state_tensor, dims=other_sites), dims=Tuple(other_sites))
    
    # Convert back to vector
    return vec(result)
end

function CMI(state::DiagonalState, A::AbstractVector, B::AbstractVector, C::AbstractVector, samples::Int; cutoff=1E-8, maxdim=200)
    S_B, S2_B = entropy(state, B, samples; cutoff=cutoff, maxdim=maxdim)
    S_AB, S2_AB = entropy(state, vcat(A,B), samples; cutoff=cutoff, maxdim=maxdim)
    S_BC, S2_BC = entropy(state, vcat(B,C), samples; cutoff=cutoff, maxdim=maxdim)
    S_ABC, S2_ABC = entropy(state, vcat(A,B,C), samples; cutoff=cutoff, maxdim=maxdim)

    CMI = S_AB + S_BC - S_B - S_ABC

    return CMI, S2_AB + S2_BC + S2_B + S2_ABC
end

function exact_CMI(state::DiagonalState, i::Int, j::Int)
    ρABC = reduced(state, i, j)
    ρAB = reduced(ρABC, 1:j-i)
    ρBC = reduced(ρABC, 2:j-i+1)
    ρB = reduced(ρABC, 2:j-i)

    # println("SAB ", entropy(ρAB), " SBC ", entropy(ρBC), " SB ", entropy(ρB), " SABC ", entropy(ρABC))

    return entropy(ρAB) + entropy(ρBC) - entropy(ρB) - entropy(ρABC)
end

function vector_CMI(ρ::AbstractVector, i::Int)
    A = 1
    B = i-1
    C = Int(log2(length(ρ))) - i
    SABC = sum([real(val) > 0.0 ? -val * log(val) : 0.0 for val in ρ])
    ρAB = sum([ρ[i:2^C:end] for i in 1:2^C])
    ρB = ρAB[1:2^B] + ρAB[2^B+1:2^(B+1)]
    ρBC = ρ[1:2^(B+C)] + ρ[2^(B+C)+1:2^(B+C+1)]

    SAB = sum([real(val) > 0.0 ? -val * log(val) : 0.0 for val in ρAB])
    SBC = sum([real(val) > 0.0 ? -val * log(val) : 0.0 for val in ρBC])
    SB = sum([real(val) > 0.0 ? -val * log(val) : 0.0 for val in ρB])

    # println("SAB ", SAB, " SBC ", SBC, " SB ", SB, " SABC ", SABC)

    return SAB + SBC - SB - SABC
end
