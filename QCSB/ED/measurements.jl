function expval(state::PureStateED, M::AbstractMatrix, pos::Int)
    ψ = state.vec
    N = Int(log2(length(ψ)))
    M_width = Int(log2(size(M)[1]))

    padded_M = pad(M, pos-1, N+1-pos-M_width)
    val = (ψ' * padded_M * ψ)/norm(ψ)^2
    return val
end

function expval(state::DiagonalStateED, M::AbstractMatrix, pos::Int)
    ρ = state.vec
    N = Int(log2(length(ρ)))
    M_width = Int(log2(size(M)[1]))

    padded_M = kron(ones(2^(pos-1)), M*ones(2^M_width), ones(2^(N+1-pos-M_width)))
    val = (padded_M' * ρ)/(sum(ρ))
    return val
end

function expval(state::MixedStateED, M::AbstractMatrix, pos::Int)
    ρ = state.mat
    N = Int(log2(size(ρ)[1]))

    M_width = Int(log2(size(M)[1]))
    padded_M = pad(M, pos-1, N+1-pos-M_width)
    val = tr(ρ * padded_M)/tr(ρ)
    return val
end

function measure(state::PureStateED, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool)
    ψ = state.vec
    N = Int(log2(length(ψ)))
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))

    padded_Π = pad(Π, pos-1, N+1-pos-M_width)
    ψ = padded_Π*ψ
    return PureStateED(ψ) / norm(PureStateED(ψ))
end

function measure(state::DiagonalStateED, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool)
    ρ = state.vec
    N = Int(log2(length(ρ)))
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    padded_Π = pad(Π, pos-1, N+1-pos-M_width)

    ρ = padded_Π*padded_Π*ρ
    return DiagonalStateED(ρ) / norm(DiagonalStateED(ρ))
end

function measure(state::MixedStateED, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool)
    ρ = state.mat
    N = Int(log2(size(ρ)[1]))
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    padded_Π = pad(Π, pos-1, N+1-pos-M_width)

    ρ =  padded_Π*ρ*padded_Π
    return MixedStateED(ρ) / norm(MixedStateED(ρ))
end