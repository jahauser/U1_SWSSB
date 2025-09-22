# warning: does not accomodate periodic boundary conditions
function expval(state::DiagonalState, M::AbstractMatrix, pos::Int; cutoff=1E-8, maxdim=200)
    ψ = state.mps
    sites = siteinds(ψ)
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[pos+i] for i in 0:M_width-1]...), ψ; cutoff=cutoff, maxdim=maxdim)
    val = inner(MPS(sites, i -> "+"), Mψ) / inner(MPS(sites, i -> "+"), ψ) 
    return val
end

function measure(state::DiagonalState, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200)
    ψ = state.mps
    sites = siteinds(ψ)
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g = op(Π, [sites[pos+i] for i in 0:M_width-1]...)

    ψ = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    return DiagonalState(ψ)
end

function measure(state::DiagonalState, M::AbstractMatrix, λ::Float64, pos::Int; cutoff=1E-8, maxdim=200)
    val = expval(state, M, pos; cutoff=cutoff, maxdim=maxdim)
    prob = (1 + 2λ/(1+λ^2)*val)/2  

    if rand() < prob
        return measure(state, M, λ, pos, false; cutoff=cutoff, maxdim=maxdim), false, val
    else
        return measure(state, M, λ, pos, true; cutoff=cutoff, maxdim=maxdim), true, val
    end
end

function forced_measure_with_val(state::DiagonalState, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200)
    return measure(state, M, λ, pos, m; cutoff=cutoff, maxdim=maxdim), m, expval(state, M, pos; cutoff=cutoff, maxdim=maxdim)
end

function measure(state::DiagonalState, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; cutoff=1E-8, maxdim=200)
    for (i,pos) in enumerate(positions)
        state = measure(state, M, λ, pos, ms[i]; cutoff=cutoff, maxdim=maxdim)
    end
    return state
end

function measure(state::DiagonalState, M::AbstractMatrix, λ::Float64, positions::AbstractVector; cutoff=1E-8, maxdim=200)
    ms = zeros(Bool, length(positions))
    vals = zeros(Float64, length(positions))

    for (i,pos) in enumerate(positions)
        state, m, val = measure(state, M, λ, pos; cutoff=cutoff, maxdim=maxdim)
        ms[i] = m
        vals[i] = val
    end
    return state, ms, vals
end

function forced_measure_with_val(state::DiagonalState, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; cutoff=1E-8, maxdim=200)
    vals = zeros(Float64, length(positions))

    for (i,pos) in enumerate(positions)
        state, _, val = forced_measure_with_val(state, M, λ, pos, ms[i]; cutoff=cutoff, maxdim=maxdim   )
        vals[i] = val
    end
    return state, ms, vals
end