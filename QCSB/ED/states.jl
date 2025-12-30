import Base: /
import ITensorMPS: truncate!, norm

# could incorporate local dimension def in here, or even system size

struct PureStateED
    vec::SparseVector{ComplexF64,Int64}
end

struct DiagonalStateED
    vec::SparseVector{ComplexF64,Int64}
end

struct MixedStateED
    mat::SparseMatrixCSC{ComplexF64,Int64}
end

StateED = Union{PureStateED, DiagonalStateED, MixedStateED}

function norm(state::PureStateED)
    return state.vec' * state.vec
end

function norm(state::DiagonalStateED)
    return sum(state.vec)
end

function norm(state::MixedStateED)
    return tr(state.mat)
end

/(state::PureStateED, scalar::Number) = PureStateED(state.vec / scalar)
/(state::DiagonalStateED, scalar::Number) = DiagonalStateED(state.vec / scalar)
/(state::MixedStateED, scalar::Number) = MixedStateED(state.mat / scalar)

truncate!(state::PureStateED; kwargs...) = state
truncate!(state::DiagonalStateED; kwargs...) = state
truncate!(state::MixedStateED; kwargs...) = state