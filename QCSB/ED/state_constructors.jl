function zero_state(::Type{PureStateED}, L::Int)
    ψ = zeros(ComplexF64, 2^L)
    ψ[1] = 1.0 + 0im
    return PureStateED(ψ)
end

function zero_state(::Type{DiagonalStateED}, L::Int)
    ρ = zeros(ComplexF64, 2^L)
    ρ[1] = 1.0 + 0im
    return DiagonalStateED(ρ)
end

function zero_state(::Type{MixedStateED}, L::Int)
    ρ = zeros(ComplexF64, 2^L, 2^L)
    ρ[1,1] = 1.0 + 0im
    return MixedStateED(ρ)
end

function ghz_state(::Type{PureStateED}, L::Int; ref=false)
    N = L+ref
    ψ = zeros(ComplexF64, 2^N)
    ψ[1] = 1.0 + 0im
    ψ[end] = 1.0 + 0im
    return PureStateED(ψ/sqrt(2))
end

function ghz_state(::Type{MixedStateED}, L::Int; ref=false)
    N = L+ref
    ρ = zeros(ComplexF64, 2^N, 2^N)
    ρ[1,1] = 1.0 + 0im
    ρ[1,end] = 1.0 + 0im
    ρ[end,1] = 1.0 + 0im
    ρ[end,end] = 1.0 + 0im
    return MixedStateED(ρ/2)
end

function block_pattern(L::Int, a::Int)
    @assert iseven(L) "L must be even"
    nblocks = Int(ceil(L / (2a)))
    v = repeat([zeros(Int, a); ones(Int, a)], nblocks)
    trim = (2a * nblocks - L) ÷ 2
    return v[trim+1:end-trim]
end

function block_state(::Type{PureStateED}, L::Int, a::Int)
    pattern = block_pattern(L, a)
    ind = parse(Int, join(pattern), base=2)
    ψ = zeros(ComplexF64, 2^L)
    ψ[ind] = 1.0 + 0im
    return PureStateED(ψ)
end

function block_state(::Type{DiagonalStateED}, L::Int, a::Int)
    pattern = block_pattern(L, a)
    ind = parse(Int, join(pattern), base=2)
    ρ = zeros(ComplexF64, 2^L)
    ρ[ind] = 1.0 + 0im
    return DiagonalStateED(ρ)
end

function block_state(::Type{MixedStateED}, L::Int, a::Int)
    pattern = repeat(block_pattern(L, a), inner=2)
    ind = parse(Int, join(pattern), base=2)
    ρ = zeros(ComplexF64, 2^L, 2^L)
    ρ[ind,ind] = 1.0 + 0im
    return MixedStateED(ρ)
end

neel_state(::Type{T}, L::Int) where T<:StateED = block_state(T, L, 1)