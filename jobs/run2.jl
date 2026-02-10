using ArgParse
using Printf
using Dates
using JLD2
using ITensors
using ITensorMPS
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "QCSB", "QCSB.jl"))

include(joinpath(ROOT, "src", "circuit.jl"))
include(joinpath(ROOT, "src", "tci.jl"))


function build_parser()
    s = ArgParseSettings(; description = "Run circuit() and save results as JLD2.")
    @add_arg_table s begin
        "--L"
            help = "System size"
            arg_type = Int
            required = true
        "--T"
            help = "Time steps"
            arg_type = Int
            required = true
        "--gamma"
            help = "Circuit parameter gamma"
            arg_type = Float64
            required = true
        "--B"
            help = "Number of sites to measure in the middle"
            arg_type = Int
            default = 0
        "--samples"
            help = "Number of samples"
            arg_type = Int
            default = 1
    end
    return s
end

# Filename-safe run tag
function tag(; L, T, gamma, B, samples)
    f(x) = replace(@sprintf("%.3f", x), "." => "p")
    return "L$(L)_T$(T)_gamma$(f(gamma))_B$(B)_samples$(samples)"
end

function simple_circuit(state::DiagonalStateMPS, L::Int, T::Int, p::Float64; cutoff=1e-8, maxdim=200)
    SWAPn1 = decoherence_layer(state, SWAP, p, 1:2:L-1)
    SWAPn2 = decoherence_layer(state, SWAP, p, 2:2:L-1)

    for t in 1:T
        state = apply(SWAPn1, state; cutoff=cutoff, maxdim=maxdim)
        state = apply(SWAPn2, state; cutoff=cutoff, maxdim=maxdim)
        state = normalize(state)
        truncate!(state; cutoff=cutoff, maxdim=maxdim)
    end

    return state
end

function number_distribution_sample(state::DiagonalStateMPS, A::Int, B::Int)
    state = deepcopy(state)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)

    state, _ = traceright(state, A+B+1)
    if B > 0
        state, _, _ = measure(state, PauliZ, 1.0, A+1:A+B; cutoff=1e-8, maxdim=200)
        state, _ = traceright(state, A+1)
    end


    state = normalize(state)

    # return state
    ψ = state.mps

    dist = ITensor(1.0)
    dist *= ψ[1]
    i1 = only(inds(ψ[1],"Site"))
    for (counter, site) in enumerate(1:A-1)
        i2 = only(inds(ψ[site+1],"Site"))
        i3 = siteind("Qudit", 1; dim=counter+2, conserve_qns=true)

        T = ITensor(dag(i1), dag(i2), i3)
        for a in 1:counter+1, b in 1:2
            c = a + b - 1
            T[i1 => a, i2 => b, i3 => c] = 1.0
        end

        dist *= ψ[site+1]*T
        i1 = i3
    end

    return Array(dist, inds(dist)...)
end

function main(args)
    opts = parse_args(args, build_parser())
    L    = opts["L"]
    T    = opts["T"]
    gamma    = opts["gamma"]
    B    = opts["B"]
    samples = opts["samples"]

    # Hardcoded output dir for easy local testing
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)
        
    _ = simple_circuit(neel_state(DiagonalStateMPS, 2; conserve_number=true), 2, 2, 0.1)
     
    t1 = time()
    state = neel_state(DiagonalStateMPS, L; conserve_number=true)
    state = simple_circuit(state, L, T, exp(-gamma)*sinh(gamma))
    
    total = 0.0
    total2 = 0.0
    for _ in 1:samples
        prob = max(abs.(number_distribution_sample(state, L÷2, B))...)
        total += prob
        total2 += prob^2
    end
    total /= samples
    total2 /= samples

    dt = time() - t1

    tagstr  = tag(L = L, T = T, gamma = gamma, B = B, samples=samples)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))

    fname = joinpath(outdir, "optimal_decoder_$(tagstr)_$(timestr)_$(randtag).jld2")

    @save fname L T gamma B samples total total2 dt
end

main(ARGS)
