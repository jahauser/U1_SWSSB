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
        "--p"
            help = "Circuit parameter p"
            arg_type = Float64
            required = true
        "--mode"
            help = "Mode: renyi2, renyi1_TCI, renyi1_sampled"
            arg_type = String
            required = true
        "--t"
            help = "Time step for renyi1 modes"
            arg_type = Int
            default = 1
        "--samples"
            help = "Number of samples for renyi1_sampled mode"
            arg_type = Int
            default = 1
    end
    return s
end

# Filename-safe run tag
function tag(; L, T, p, mode, t, samples)
    f(x) = replace(@sprintf("%.3f", x), "." => "p")
    return "L$(L)_T$(T)_p$(f(p))_mode$(mode)_t$(t)_samples$(samples)"
end

function main(args)
    opts = parse_args(args, build_parser())
    L    = opts["L"]
    T    = opts["T"]
    p    = opts["p"]
    mode = opts["mode"]
    t    = opts["t"]
    samples = opts["samples"]

    # Hardcoded output dir for easy local testing
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)

    if mode == "renyi2"
        observables = [:κ2, :inner_κ2, :κ2_corrs]
        data_ts = vcat([0], 1:10:T)
    elseif mode == "renyi1_TCI"
        observables = [:κ1_TCI, :inner_κ1_TCI, :κ1_corrs_TCI]
        data_ts = [t]
    elseif mode == "renyi1_sampled"
        observables = [:κ1_corrs_sampled]
        data_ts = [t]
    end
        
    _, _ = circuit(neel_state(DiagonalStateMPS, 2), 2, 2, 0.1; observables=observables, renyi1_samples=samples, data_ts=data_ts)
     
    t1 = time()
    state = neel_state(DiagonalStateMPS, L)
    state, data = circuit(state, L, T, p; observables=observables, renyi1_samples=samples, data_ts=data_ts)
    dt = time() - t1

    tagstr  = tag(L = L, T = T, p = p, mode=mode, t=t, samples=samples)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))

    fname = joinpath(outdir, "sample_$(tagstr)_$(timestr)_$(randtag).jld2")

    @save fname L T p mode t samples observables data dt
end

main(ARGS)
