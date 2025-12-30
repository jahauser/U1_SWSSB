function measure(state::State, M::AbstractMatrix, λ::Float64, pos::Int; kwargs...)
    val = expval(state, M, pos; kwargs...)
    prob = (1 + 2λ/(1+λ^2)*val)/2  

    if rand() < abs(prob)
        return measure(state, M, λ, pos, false; kwargs...), false, val
    else
        return measure(state, M, λ, pos, true; kwargs...), true, val
    end
end

function forced_measure_with_val(state::State, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; kwargs...)
    return measure(state, M, λ, pos, m; kwargs...), m, expval(state, M, pos; kwargs...)
end

function measure(state::State, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; kwargs...)
    for (i,pos) in enumerate(positions)
        state = measure(state, M, λ, pos, ms[i]; kwargs...)
    end
    return state
end

function measure(state::State, M::AbstractMatrix, λ::Float64, positions::AbstractVector; kwargs...)
    ms = zeros(Bool, length(positions))
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        state, m, val = measure(state, M, λ, pos; kwargs...)
        ms[i] = m
        vals[i] = val
    end
    return state, ms, vals
end

function forced_measure_with_val(state::State, M::AbstractMatrix, λ::Float64, positions::AbstractVector, ms::Vector{Bool}; kwargs...)
    vals = zeros(ComplexF64, length(positions))

    for (i,pos) in enumerate(positions)
        state, _, val = forced_measure_with_val(state, M, λ, pos, ms[i]; kwargs...)
        vals[i] = val
    end
    return state, ms, vals
end