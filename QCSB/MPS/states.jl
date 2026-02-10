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

# function norm(state::DiagonalStateMPS)
#     ψ = state.mps
#     sites = siteinds(ψ)
#     L = length(sites)
#     return exp(loginner(MPS(sites, i -> "+"), ψ))*2^(L/2)
# end

# function 

function traceright(state::DiagonalStateMPS, R::Int)
    ψ = copy(state.mps)
    sites = siteinds(ψ)
    L = length(sites)

    if R > L
        return state, 0.0
    end

    mats = Vector{ITensor}(undef, L-R+1)
    for j in R:L
        s = sites[j]
        mats[j-R+1] = dag(ITensors.state(s, "0") + ITensors.state(s, "1")) * ψ[j]
    end

    v₀ = mats[end]
    logscale = 0.0
    for j in 1:(L-R)
        v₁ = mats[end-j] * v₀
        n = norm(v₁)
        logscale += log(n)
        v₀ = v₁ / n
    end

    if R == 1
        return nothing, logscale + log(Complex(Array(v₀)[]))
    else
    
        ψ[R-1] *= v₀

        s = exp(+ logscale / (L - R + 1))
        for j in 1:(L-R+1)
            ψ[j] *= s
        end
    return DiagonalStateMPS(MPS(ψ[1:R-1])), logscale
    end
end

function normalize(state::DiagonalStateMPS)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    _, λ = traceright(state, 1)
    λ /= L
    
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