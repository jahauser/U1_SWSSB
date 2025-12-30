import Base: /
import ITensorMPS: truncate!, norm, apply

struct PureStateMPS
    mps::MPS
end

struct DiagonalStateMPS
    mps::MPS
end

struct MixedStateMPS
    mps::MPS
end

StateMPS = Union{PureStateMPS, DiagonalStateMPS, MixedStateMPS}


apply(gates::Vector, state::T; kwargs...) where {T <: StateMPS} = T(apply(gates, state.mps; kwargs...))

function norm(state::PureStateMPS)
    ψ = state.mps
    return norm(ψ)
end

# function norm(state::DiagonalStateMPS)
#     ψ = state.mps
#     sites = siteinds(ψ)
#     L = length(sites)
#     return inner(MPS(sites, i -> "+"), ψ)*2^(L/2)
# end

function norm(state::DiagonalStateMPS)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    return exp(loginner(MPS(sites, i -> "+"), ψ))*2^(L/2)
end

function normalize(state::DiagonalStateMPS)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    λ = (loginner(MPS(sites, i -> "+"), ψ) + (L/2)*log(2)) / L
    s = exp(-λ)
    for j in 1:L
        ψ[j] *= s
    end
    return DiagonalStateMPS(ψ)
end

function norm(state::MixedStateMPS)
    ψ = state.mps
    return doubledtrace(ψ)
end

/(state::PureStateMPS, scalar::Number) = PureStateMPS(state.mps / scalar)
truncate!(state::PureStateMPS; kwargs...) = PureStateMPS(truncate!(state.mps; kwargs...))

/(state::DiagonalStateMPS, scalar::Number) = DiagonalStateMPS(state.mps / scalar)
truncate!(state::DiagonalStateMPS; kwargs...) = DiagonalStateMPS(truncate!(state.mps; kwargs...))

/(state::MixedStateMPS, scalar::Number) = MixedStateMPS(state.mps / scalar)
truncate!(state::MixedStateMPS; kwargs...) = MixedStateMPS(truncate!(state.mps; kwargs...))

dense(state::PureStateMPS) = to_vector(state.mps)
dense(state::DiagonalStateMPS) = spdiagm(to_vector(state.mps))
dense(state::MixedStateMPS) = to_matrix(state.mps)

function to_vector(ψ::MPS)
    return to_vector(contract(ψ))
end

function to_vector(ψ::ITensor)
    tensor_inds = inds(ψ)
    permutation = vcat(reverse(tensor_inds))

    M = Array{ComplexF64,length(tensor_inds)}(ψ, permutation)

    return reshape(M, 2^(length(tensor_inds)))
end

function to_matrix(ψ::MPS)
    return to_matrix(contract(ψ))
end

function to_matrix(ψ::ITensor)
    tensor_inds = inds(ψ)
    permutation = vcat(reverse(tensor_inds[2:2:end]), reverse(tensor_inds[1:2:end]))

    M = Array{ComplexF64,length(tensor_inds)}(ψ, permutation)

    return reshape(M, (2^(length(tensor_inds)÷2), 2^(length(tensor_inds)÷2)))
end