# (1)0101...01 = (4^(n+1)-1)/3
function MPS_neelstate(L::Int)
    N = L
    sites = siteinds("Qubit", N)
    ρ = MPS(sites, i -> mod(i-1+L,2)==0 ? "0" : "1")
    return ρ, sites
end

# function MPS_dwstate(L::Int)
#     N = L
#     sites = siteinds("Qubit", N)
#     ρ = MPS(sites, i -> i <= L÷2 ? "0" : "1")
#     return ρ, sites
# end

function to_vector(ρ::MPS)
    ρ = contract(ρ)
    tensor_inds = inds(ρ)
    permutation = vcat(reverse(tensor_inds))

    M = Array{ComplexF64,length(tensor_inds)}(ρ, permutation)

    return reshape(M, 2^(length(tensor_inds)))
end

function to_vector(ρ::ITensor)
    tensor_inds = inds(ρ)
    permutation = vcat(reverse(tensor_inds))

    M = Array{ComplexF64,length(tensor_inds)}(ρ, permutation)

    return reshape(M, 2^(length(tensor_inds)))
end





import Base: /
import ITensorMPS: truncate!

struct DiagonalState
    mps::MPS
end

function diagonal_neelstate(L::Int)
    mps = MPS_neelstate(L)[1]
    state = DiagonalState(mps)
    return state / norm(state)
end

function norm(state::DiagonalState)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    return inner(MPS(sites, i -> "+"), ψ)*2^(L/2)
end

/(state::DiagonalState, scalar::Number) = DiagonalState(state.mps / scalar)

truncate!(state::DiagonalState; kwargs...) = DiagonalState(truncate!(state.mps; kwargs...))