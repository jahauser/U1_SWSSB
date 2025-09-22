function imaginary_time_evolve(state::DiagonalState, h::AbstractMatrix, γ::Float64, positions::AbstractVector; cutoff=1E-8, maxdim=200)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)
    h_width = Int(log2(size(h)[1]))

    gates = ITensor[]
    for pos in positions
        M = op(exp(γ * h), [sites[mod1(pos+i,L)] for i in 0:h_width-1]...)
        push!(gates, M)
    end

    ψ = apply(gates, ψ; cutoff=cutoff, maxdim=maxdim)
    
    # Normalize the state
    ψ /= inner(MPS(sites, i -> "+"), ψ)*2^(length(ψ)/2)
    
    # Truncate the MPS
    truncate!(ψ; cutoff=cutoff, maxdim=maxdim)

    return DiagonalState(ψ)
end

# function imaginary_time_evolve(ψ::MPS, h::AbstractMatrix, γ::Float64, positions::AbstractVector; cutoff=1E-8, maxdim=200)
#     sites = siteinds(ψ)
#     h_width = Int(log2(size(h)[1]))

#     gates = ITensor[]
#     for pos in positions
#         M = op(exp(γ * h), [sites[pos+i] for i in 0:h_width-1]...)
#         push!(gates, M)
#     end

#     ψ = apply(gates, ψ; cutoff=cutoff, maxdim=maxdim)
    
#     # Normalize the state
#     ψ /= inner(MPS(sites, i -> "+"), ψ)*2^(length(ψ)/2)
    
#     # Truncate the MPS
#     truncate!(ψ; cutoff=cutoff, maxdim=maxdim)
    
#     return ψ
# end