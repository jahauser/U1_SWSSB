function zero_state(::Type{PureStateMPS}, L::Int)
    sites = siteinds("Qubit", L)
    ψ = MPS(sites, _ -> "0")
    return PureStateMPS(ψ)
end

function zero_state(::Type{DiagonalStateMPS}, L::Int)
    sites = siteinds("Qubit", L)
    ψ = MPS(sites, _ -> "0")
    return DiagonalStateMPS(ψ)
end

function zero_state(::Type{MixedStateMPS}, L::Int)
    sites = siteinds("Qubit", 2L)
    ψ = MPS(sites, _ -> "0")
    return MixedStateMPS(ψ)
end

function ghz_state(::Type{PureStateMPS}, L::Int; ref=false)
    N = L+ref
    sites = siteinds("Qubit", N)
    ψ0 = MPS(sites, _ -> "0")
    ψ1 = MPS(sites, _ -> "1")
    return PureStateMPS((ψ0 + ψ1)/sqrt(2))
end

function ghz_state(::Type{MixedStateMPS}, L::Int; ref=false)
    N = 2L+2ref
    sites = siteinds("Qubit", N)
    ρ00 = MPS(sites, _ -> "0")
    ρ11 = MPS(sites, _ -> "1")
    ρ01 = MPS(sites, i -> mod(i,2)==0 ? "0" : "1")
    ρ10 = MPS(sites, i -> mod(i,2)==0 ? "1" : "0")
    return MixedStateMPS((ρ00 + ρ01 + ρ10 + ρ11)/2)
end

function block_pattern(L::Int, a::Int)
    @assert iseven(L) "L must be even"
    nblocks = Int(ceil(L / (2a)))
    v = repeat([zeros(Int, a); ones(Int, a)], nblocks)
    trim = (2a * nblocks - L) ÷ 2
    return v[trim+1:end-trim]
end

function block_state(::Type{PureStateMPS}, L::Int, a::Int; conserve_number=false)
    pattern = block_pattern(L, a)
    sites = siteinds("Qubit", L; conserve_number=conserve_number)
    ψ = MPS(sites, i -> pattern[i] == 0 ? "0" : "1")
    return PureStateMPS(ψ)
end

function block_state(::Type{DiagonalStateMPS}, L::Int, a::Int; conserve_number=false)
    pattern = block_pattern(L, a)
    sites = siteinds("Qubit", L; conserve_number=conserve_number)
    ψ = MPS(sites, i -> pattern[i] == 0 ? "0" : "1")
    return DiagonalStateMPS(ψ)
end

function block_state(::Type{MixedStateMPS}, L::Int, a::Int; conserve_number=false)
    pattern = repeat(block_pattern(L, a), inner=2)
    sites = siteinds("Qubit", 2L; conserve_number=conserve_number)
    ψ = MPS(sites, i -> pattern[i] == 0 ? "0" : "1")
    return MixedStateMPS(ψ)
end

neel_state(::Type{T}, L::Int; conserve_number=false) where T<:StateMPS = block_state(T, L, 1; conserve_number=conserve_number)