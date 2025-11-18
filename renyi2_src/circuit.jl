Id = [1 0 
      0 1]

PauliX = [0 1
          1 0]

PauliY = [0 -1im
          1im 0]

PauliZ = [1 0
          0 -1]

SWAP = [1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1]

ITensors.op(::OpName"a",::SiteType"Qubit") = [0 1
                                              0 0]
ITensors.op(::OpName"a†",::SiteType"Qubit") = [0 0
                                              1 0]

function χ2_corr(state::DiagonalStateMPS)
    state /= norm(state)
    return correlation_matrix(state.mps, "a", "a†")
end

function χ2(state::DiagonalStateMPS)
    L = length(state.mps)

    corrs = χ2_corr(state)
    corrs[diagind(corrs)] .= 0
    χ = sum(corrs)/L
    return χ
end

function initiate_data(state::DiagonalStateMPS, observables::Vector{Symbol}, L::Int, T::Int)
    data = Dict{Symbol,Array}()
    if :χ2 in observables
        data[:χ2] = zeros(Float64, T+1)
    end
    return data
end

function update_data(state::DiagonalStateMPS, observables::Vector{Symbol}, data::Dict{Symbol,Array}, t::Int)
    if :χ2 in observables
        data[:χ2][t] = χ2(state)
    end
    return data
end

function circuit(state::DiagonalStateMPS, L::Int, T::Int, p::Float64; observables=Symbol[], cutoff=1E-8, maxdim=200)
    SWAPn1 = decoherence_layer(state, SWAP, p, 1:2:L-1)
    SWAPn2 = decoherence_layer(state, SWAP, p, 2:2:L-1)

    data = initiate_data(state, observables, L, T)
    data = update_data(state, observables, data, 1)

    for t in 1:T
        state = apply(SWAPn1, state; cutoff=cutoff, maxdim=maxdim)
        state = apply(SWAPn2, state; cutoff=cutoff, maxdim=maxdim)
        state /= norm(state)
        truncate!(state; cutoff=cutoff, maxdim=maxdim)

        data = update_data(state, observables, data, t+1)
    end

    return state, data
end