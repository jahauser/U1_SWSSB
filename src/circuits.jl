function initiate_data(_::DiagonalState, observables::Vector{Symbol}, L::Int, T::Int)
    data = Dict{Symbol,Vector}()
    # if :κ2_ratio in observables
    #     data[:κ2_ratio] = zeros(Float64, 2T+2)
    # end
    # if :maxlinkdim in observables
    #     data[:maxlinkdim] = zeros(Float64, 2T+2)
    # end
    # if :S in observables
    #     data[:S] = zeros(Float64, 2T+2)
    # end
    # if :Z_profile in observables
    #     data[:Z_profile] = [zeros(Float64, L) for _ in 1:2T+2]
    # end
    # if :Z2_profile in observables
    #     data[:Z2_profile] = [zeros(Float64, L, L) for _ in 1:2T+2]
    # end
    # if :CMI in observables
    #     data[:CMI] = [zeros(Float64, L-2) for _ in 1:2T+2]
    # end
    # if :mid_CMI in observables
    #     data[:mid_CMI] = zeros(Float64, 2T+2)
    # end
    # if :exact_CMI in observables
    #     data[:exact_CMI] = [zeros(Float64, L-2) for _ in 1:2T+2]
    # end
    # if :exact_mid_CMI in observables
    #     data[:mid_CMI] = zeros(Float64, 2T+2)
    # end
    if :exact_CMI_corr in observables
        data[:exact_CMI_corr] = zeros(Float64, 2T+2)
    end
    return data
end

function update_data(state::DiagonalState, observables::Vector{Symbol}, data::Dict{Symbol,Vector}, t::Int; r=2)
    ψ = state.mps
    L = length(ψ)
    # if :κ2_ratio in observables
    #     corrs = χ2_corr(ψ)
    #     i1 = L÷4
    #     j1 = 3L÷4+1
    #     i2 = 3L÷8
    #     j2 = 5L÷8+1
    #     data[:κ2_ratio][t] = corrs[i1,j1]/corrs[i2,j2]
    # end
    # if :maxlinkdim in observables
    #     data[:maxlinkdim][t] = maxlinkdim(ψ)
    # end
    # if :S in observables
    #     data[:S][t] = entanglemententropy(ψ)
    # end
    # if :Z_profile in observables
    #     data[:Z_profile][t] = 2expect(ψ, "Sz")
    # end
    # if :Z2_profile in observables
    #     Zs = 2expect(ψ, "Sz")
    #     data[:Z2_profile][t] = Zs * Zs'
    # end
    # if :CMI in observables
    #     samples = 1000
    #     CMIs = zeros(L-2)
    #     for i in 2:L-1
    #         CMIs[i-1], _ = CMI(state, 1:1, 2:i, i+1:L, samples; cutoff=1E-8, maxdim=200)
    #     end
    #     data[:CMI][t] = CMIs
    # end
    # if :mid_CMI in observables
    #     samples = 1000
    #     data[:mid_CMI][t], _ = CMI(state, 1:1, 2:L÷2, L÷2+1:L, samples; cutoff=1E-8, maxdim=200)
    # end
    # if :exact_CMI in observables
    #     ρ = to_vector(ψ)
    #     CMIs = zeros(L-2)
    #     for i in 2:L-1
    #         CMIs[i-1] = vector_CMI(ρ, i)
    #     end
    #     data[:exact_CMI][t] = CMIs
    # end
    # if :exact_mid_CMI in observables
    #     ρ = to_vector(ψ)
    #     data[:exact_mid_CMI][t] = vector_CMI(ρ, L÷2)
    # end
    if :exact_CMI_corr in observables
        data[:exact_CMI_corr][t] = exact_CMI(state, L÷2+1-r, L÷2+1+r)
    end
    return data
end


function circuit(state::DiagonalState, L::Int, T::Int, γ::Float64, λ::Float64; observables=Symbol[], terminal_observables=Symbol[],
     PBC=false, cutoff=1E-8, maxdim=200, r=r)

    data = initiate_data(state, observables, L, T)
    data = update_data(state, observables, data, 1; r=r)
    data = update_data(state, observables, data, 2; r=r)

    terminal_data = initiate_data(state, terminal_observables, L, 2)

    for t in 1:T
        state = imaginary_time_evolve(state, SWAP, γ, 1:2:L-1+PBC; cutoff=cutoff, maxdim=maxdim)
        state = imaginary_time_evolve(state, SWAP, γ, 2:2:L-1+PBC; cutoff=cutoff, maxdim=maxdim)
        state /= norm(state)
        truncate!(state; cutoff=cutoff, maxdim=maxdim)
        data = update_data(state, observables, data, 2t+1; r=r)

        if t == T
            terminal_data = update_data(state, terminal_observables, terminal_data, 1; r=r)
        end

        if λ > 0.0
            state, _, _ = measure(state, PauliZ, λ, 1:L; cutoff=cutoff, maxdim=maxdim)
            state /= norm(state)
            truncate!(state; cutoff=cutoff, maxdim=maxdim)
            data = update_data(state, observables, data, 2t+2; r=r)
        end

        if t == T
            terminal_data = update_data(state, terminal_observables, terminal_data, 2; r=r)
        end
    end

    return state, data, terminal_data
end

function circuit(state::DiagonalState, L::Int, T::Int, γ::Float64, λ::Float64, samples::Int; observables=Symbol[],
    terminal_observables=Symbol[], PBC=false, cutoff=1E-8, maxdim=200, r=r)

    full_data = initiate_data(state, observables, L, T)
    full_terminal_data = initiate_data(state, terminal_observables, L, 2)

    for _ in 1:samples
        _, data, terminal_data = circuit(state, L, T, γ, λ; observables=observables, PBC=PBC, cutoff=cutoff, maxdim=maxdim, r=r)
        for observable in observables
            full_data[observable] .+= data[observable]
        end
        for observable in terminal_observables
            full_terminal_data[observable] .+= terminal_data[observable]
        end
    end

    for observable in observables
        full_data[observable] /= samples
    end
    for observable in terminal_observables
        full_terminal_data[observable] /= samples
    end
    return full_data, full_terminal_data
end


# function MPS_initiate_data(observables::Vector{Symbol}, L::Int, T::Int)
#     data = Dict{Symbol,Vector}()
#     if :χ2 in observables
#         data[:χ2] = zeros(Float64, T+1)
#     end
#     if :χ2_corr in observables
#         data[:χ2_corr] = zeros(Float64, T+1)
#     end
#     if :maxlinkdim in observables
#         data[:maxlinkdim] = zeros(Float64, T+1)
#     end
#     if :S in observables
#         data[:S] = zeros(Float64, T+1)
#     end
#     if :Z_profile in observables
#         data[:Z_profile] = [zeros(Float64, L) for _ in 1:T+1]
#     end
#     return data
# end

# function MPS_update_data(ψ::MPS, observables::Vector{Symbol}, data::Dict{Symbol,Vector}, t::Int)
#     L = length(ψ)
#     if :χ2 in observables
#         data[:χ2][t] = χ2(ψ)
#     end
#     if :χ2_corr in observables
#         data[:χ2_corr][t] = χ2_corr(ψ)[L÷4,3L÷4+1]
#     end
#     if :maxlinkdim in observables
#         data[:maxlinkdim][t] = maxlinkdim(ψ)
#     end
#     if :S in observables
#         data[:S][t] = entanglemententropy(ψ)
#     end
#     if :Z_profile in observables
#         data[:Z_profile][t] = 2expect(ψ, "Sz")
#     end
#     return data
# end

# function MPS_classical(ψ::MPS, L::Int, T::Int, γ::Float64; observables=Symbol[], cutoff=1E-8, maxdim=200)
#     sites = siteinds(ψ)

#     data = MPS_initiate_data(observables, L, T)
#     data = MPS_update_data(ψ, observables, data, 1)

#     for t in 1:T
#         ψ = imaginary_time_evolve(ψ, SWAP, γ, 1:2:L-1; cutoff=cutoff, maxdim=maxdim)
#         ψ = imaginary_time_evolve(ψ, SWAP, γ, 2:2:L-1; cutoff=cutoff, maxdim=maxdim)
#         ψ /= inner(MPS(sites, i -> "+"), ψ)*2^(L/2)
#         truncate!(ψ; cutoff=cutoff, maxdim=maxdim)

#         data = MPS_update_data(ψ, observables, data, t+1)
#     end

#     return ψ, data
# end
