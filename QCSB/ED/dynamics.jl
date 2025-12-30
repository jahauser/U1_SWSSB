struct UnitaryLayer
    gates::Vector{Tuple{Int, AbstractMatrix}}
end

struct DecoherenceLayer
    gates::Vector{Tuple{Int, AbstractMatrix, Float64}}
end

function unitary_layer(_::StateED, M::AbstractMatrix, θ::Float64, positions::AbstractVector)
    return UnitaryLayer( [ (pos, exp(1im*θ*M)) for pos in positions ] )
end

function decoherence_layer(_::StateED,  M::AbstractMatrix, p::Float64, positions::AbstractVector)
    return DecoherenceLayer( [ (pos, M, p) for pos in positions ] )
end

function apply(layer::UnitaryLayer, state::PureStateED)
    ψ = state.vec
    N = Int(log2(length(ψ)))
    for (pos, gate) in layer.gates
        gate_width = Int(log2(size(gate)[1]))
        padded_M = pad(gate, pos-1, N+1-pos-gate_width)
        ψ = padded_M * ψ
    end
    return PureStateED(ψ)
end

function apply(layer::UnitaryLayer, state::MixedStateED)
    ρ = state.mat
    N = Int(log2(size(ρ)[1]))
    for (pos, gate) in layer.gates
        gate_width = Int(log2(size(gate)[1]))
        padded_M = pad(gate, pos-1, N+1-pos-gate_width)
        ρ = padded_M * ρ * padded_M
    end
    return MixedStateED(ρ)
end

# assumes M satisfies M^2 = 1
function apply(layer::DecoherenceLayer, state::DiagonalStateED)
    ρ = state.vec
    N = Int(log2(length(ρ)))
    for (pos, gate, p) in layer.gates
        gate_width = Int(log2(size(gate)[1]))
        padded_M = pad(gate, pos-1, N+1-pos-gate_width)
        ρ = (1-p)*ρ + p*padded_M*ρ
    end
    return DiagonalStateED(ρ)
end

function apply(layer::DecoherenceLayer, state::MixedStateED)
    ρ = state.mat
    N = Int(log2(size(ρ)[1]))
    for (pos, gate, p) in layer.gates
        gate_width = Int(log2(size(gate)[1]))
        padded_M = pad(gate, pos-1, N+1-pos-gate_width)
        ρ = (1-p)*ρ + p*padded_M*ρ*padded_M
    end
    return MixedStateED(ρ)
end

# function decohere(state::DiagonalStateED, M::AbstractMatrix, p::Float64, site::Int)
#     L = Int(log2(size(state.mat)[1]))
#     M_width = Int(log2(size(M)[1]))
#     padded_M = pad(M, site-1, L+1-site-M_width)
#     return DiagonalStateED( (1-p)*state.mat + p*padded_M*state.mat*padded_M )
# end

# function decohere(state::DiagonalStateED, M::AbstractMatrix, p::Float64, sites::AbstractVector)
#     for site in sites
#         state = decohere(state, M, p, site)
#     end
#     return state
# end

# # assumes M satisfies M^2 = 1
# function decohere(state::MixedStateED, M::AbstractMatrix, p::Float64, site::Int)
#     L = Int(log2(size(state.mat)[1]))
#     M_width = Int(log2(size(M)[1]))
#     padded_M = pad(M, site-1, L+1-site-M_width)
#     return MixedStateED( (1-p)*state.mat + p*padded_M*state.mat*padded_M )
# end

# function decohere(state::MixedStateED, M::AbstractMatrix, p::Float64, sites::AbstractVector)
#     for site in sites
#         state = decohere(state, M, p, site)
#     end
#     return state
# end
