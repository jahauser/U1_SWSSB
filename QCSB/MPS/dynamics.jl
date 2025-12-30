function decoherence_layer(state::DiagonalStateMPS,  M::AbstractMatrix, p::Float64, positions::AbstractVector; refs=0)
    sites = siteinds(state.mps)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        legs = [sites[mod1(pos+i,L)] for i in 0:M_width-1]
        g = (1-p)*op(I, legs...) + p*op(M, legs...)
        push!(gates, g)
    end

    return gates
end

function decoherence_layer(state::MixedStateMPS, M::AbstractMatrix, p::Float64, positions::AbstractVector; refs=0)
    sites = siteinds(state.mps)
    L = length(sites)÷2 - refs
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        legs1 = [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]
        legs2 = [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]
        g = (1-p)*op(I, legs1...)*op(I, legs2...) + p*op(M, legs1...)*op(M, legs2...)
        push!(gates, g)
    end
    return gates
end

function unitary_layer(state::PureStateMPS, M::AbstractMatrix, θ::Float64, positions::AbstractVector; refs=0)
    sites = siteinds(state.mps)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        legs = [sites[mod1(pos+i,L)] for i in 0:M_width-1]
        h = op(exp(1im*θ*M), legs...)
        push!(gates, h)
    end
    return gates
end

function unitary_layer(state::MixedStateMPS, M::AbstractMatrix, θ::Float64, positions::AbstractVector; refs=0)
    sites = siteinds(state.mps)
    L = length(sites)÷2 - refs
    M_width = Int(log2(size(M)[1]))

    gates = ITensor[]
    for pos in positions
        legs1 = [sites[mod1(2*(pos+i)-1,2L)] for i in 0:M_width-1]
        legs2 = [sites[mod1(2*(pos+i),2L)] for i in 0:M_width-1]
        h1 = op(exp(1im*θ*M), legs1...)
        h2 = op(exp(-1im*θ*M), legs2...)
        push!(gates, h1)
        push!(gates, h2)
    end
    return gates
end

