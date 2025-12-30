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

function susceptibility(corrs::AbstractMatrix)
    L = size(corrs, 1)

    corrs[diagind(corrs)] .= 0
    χ = sum(corrs)/L
    return χ
end

function inner_susceptibility(corrs::AbstractMatrix)
    L = size(corrs, 1)

    corrs[diagind(corrs)] .= 0
    χ = sum(corrs[L÷4+1:3L÷4,L÷4+1:3L÷4])/(L/2)
    return χ
end

function χ1_sweep_sample(state::DiagonalStateMPS; samples=1)
    state = normalize(state)
    L = length(state.mps)
    corrs = zeros(Float64, L÷2)
    for i in 1:L÷2
        j = L - i + 1
        sum_corr = 0.0
        for _ in 1:samples
            state_copy = deepcopy(state)
            state_copy, _, vals1 = forced_measure_with_val(state_copy, PauliZ, 1.0, [i, j], [true, false])
            state_copy, ms, vals2 = measure(state_copy, PauliZ, 1.0, [k for k in 1:L if k != i && k != j])
            if isnan(vals1[1]) || isnan(vals1[2]) || any(isnan.(vals2))
                continue
            end
            logprob1 = log((1-vals1[1])/2) + log((1+vals1[2])/2)
            logprob2 = sum([log((1 + (-1)^(ms[k]) * vals2[k])/2) for k in 1:length(ms)])
            total_logprob = logprob1 + logprob2
            total_prob1 = exp(total_logprob)


            ms2 = zeros(Bool, L)
            ms2[1:i-1] = ms[1:i-1]
            ms2[i] = false
            ms2[i+1:j-1] = ms[i:j-2]
            ms2[j] = true
            ms2[j+1:L] = ms[j-1:end]
            total_prob2 = get_prob_vec(state.mps, ms2)
            # total_prob2 = 1.0

            sum_corr += exp(logprob1) * sqrt(total_prob2 / total_prob1)
        end
        corrs[i] = abs(sum_corr / samples)
    end
    return corrs
end

function initiate_data(state::DiagonalStateMPS, observables::Vector{Symbol}, L::Int, T::Int)
    data = Dict{Symbol,Array}()
    if :κ2 in observables
        data[:κ2] = zeros(Float64, T)
    end
    if :inner_κ2 in observables
        data[:inner_κ2] = zeros(Float64, T)
    end
    if :κ2_corrs in observables
        data[:κ2_corrs] = zeros(Float64, T, L÷2)
    end
    if :κ1_TCI in observables
        data[:κ1_TCI] = zeros(Float64, T)
    end
    if :inner_κ1_TCI in observables
        data[:inner_κ1_TCI] = zeros(Float64, T)
    end
    if :κ1_corrs_TCI in observables
        data[:κ1_corrs_TCI] = zeros(Float64, T, L÷2)
    end
    if :κ1_corrs_sampled in observables
        data[:κ1_corrs_sampled] = zeros(Float64, T, L÷2)
    end
    return data
end

function update_data(state::DiagonalStateMPS, observables::Vector{Symbol}, data::Dict{Symbol,Array}, t::Int; renyi1_samples=1)
    if :κ2 in observables || :inner_κ2 in observables || :κ2_corrs in observables
        renyi2_mps = state.mps
        renyi2_mps /= exp(lognorm(renyi2_mps))
        renyi2_corrs = correlation_matrix(renyi2_mps, "a†", "a")
    end
    if :κ1_TCI in observables || :inner_κ1_TCI in observables || :κ1_corrs_TCI in observables
        renyi1_mps = tci(state)
        renyi1_mps /= exp(lognorm(renyi1_mps))
        renyi1_corrs = correlation_matrix(renyi1_mps, "a†", "a")
    end
    if :κ2 in observables
        data[:κ2][t] = susceptibility(renyi2_corrs)
    end
    if :inner_κ2 in observables
        data[:inner_κ2][t] = inner_susceptibility(renyi2_corrs)
    end
    if :κ2_corrs in observables
        L = length(renyi2_mps)
        data[:κ2_corrs][t,:] = [renyi2_corrs[i,L-i+1] for i in 1:L÷2]
    end
    if :κ1_TCI in observables
        data[:κ1_TCI][t] = susceptibility(renyi1_corrs)
    end
    if :inner_κ1_TCI in observables
        data[:inner_κ1_TCI][t] = inner_susceptibility(renyi1_corrs)
    end
    if :κ1_corrs_TCI in observables
        L = length(renyi1_mps)
        data[:κ1_corrs_TCI][t,:] = [renyi1_corrs[i,L-i+1] for i in 1:L÷2]
    end
    if :κ1_corrs_sampled in observables
        data[:κ1_corrs_sampled][t,:] = χ1_sweep_sample(state; samples=renyi1_samples)
    end
    return data
end

function circuit(state::DiagonalStateMPS, L::Int, T::Int, p::Float64; observables=Symbol[], cutoff=1E-8, maxdim=200,
                 renyi1_samples=1, data_ts=0:T)
    SWAPn1 = decoherence_layer(state, SWAP, p, 1:2:L-1)
    SWAPn2 = decoherence_layer(state, SWAP, p, 2:2:L-1)

    data = initiate_data(state, observables, L, length(data_ts))

    if 0 in data_ts
        data = update_data(state, observables, data, 1; renyi1_samples=renyi1_samples)
        data_counter = 2
    else 
        data_counter = 1
    end


    for t in 1:T
        state = apply(SWAPn1, state; cutoff=cutoff, maxdim=maxdim)
        state = apply(SWAPn2, state; cutoff=cutoff, maxdim=maxdim)
        state = normalize(state)
        truncate!(state; cutoff=cutoff, maxdim=maxdim)

        if t in data_ts
            data = update_data(state, observables, data, data_counter; renyi1_samples=renyi1_samples)
            data_counter += 1
        end

        if data_counter > length(data_ts)
            break
        end
    end

    return state, data
end