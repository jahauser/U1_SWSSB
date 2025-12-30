function expval(state::PureStateMPS, M::AbstractMatrix, pos::Int; cutoff=1E-8, maxdim=200, refs=0)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ; cutoff=cutoff, maxdim=maxdim)
    val = inner(ψ, Mψ) / norm(ψ)^2
    return val
end

function expval(state::DiagonalStateMPS, M::AbstractMatrix, pos::Int; cutoff=1E-8, maxdim=200, refs=0)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ; cutoff=cutoff, maxdim=maxdim)
    # val = inner(MPS(sites, i -> "+"), Mψ) / inner(MPS(sites, i -> "+"), ψ) 
    val = exp(loginner(MPS(sites, i -> "+"), Mψ) - loginner(MPS(sites, i -> "+"), ψ))
    return val
end

function expval(state::MixedStateMPS, M::AbstractMatrix, pos::Int; cutoff=1E-8, maxdim=200, refs=0)
    ρ = state.mps
    sites = siteinds(ρ)

    bra = M_bra(sites, M, pos; refs=refs)
    val = inner(bra, ρ) / doubledtrace(ρ)
    return val
end

function measure(state::PureStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200, refs=0)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g = op(Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    return PureStateMPS(ψ) / norm(PureStateMPS(ψ))
end

function measure(state::DiagonalStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200, refs=0)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g = op(Π*Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    return DiagonalStateMPS(ψ) / norm(DiagonalStateMPS(ψ))
end

function measure(state::MixedStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200, refs=0)
    ρ = state.mps
    sites = siteinds(ρ)
    M_width = Int(log2(size(M)[1]))
    L = length(sites)÷2 - refs

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g1 = op(Π, [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]...)
    g2 = op(Π, [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]...)

    ρ = apply([g1, g2], ρ; cutoff=cutoff, maxdim=maxdim)
    return MixedStateMPS(ρ) / norm(MixedStateMPS(ρ))
end