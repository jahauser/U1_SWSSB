using ArgParse
using Printf
using Dates
using JLD2
using ITensors
using ITensorMPS
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "qcsb", "MPS", "states.jl"))
include(joinpath(ROOT, "qcsb", "MPS", "state_constructors.jl"))
include(joinpath(ROOT, "qcsb", "MPS", "dynamics.jl"))
include(joinpath(ROOT, "qcsb", "MPS", "tools.jl"))
include(joinpath(ROOT, "qcsb", "MPS", "measurements.jl"))
include(joinpath(ROOT, "qcsb", "MPS", "entropies.jl"))

State = StateMPS
include(joinpath(ROOT, "qcsb", "src", "global", "states.jl"))

include(joinpath(ROOT, "renyi2_src", "circuit.jl"))


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
    end
    return s
end

# Filename-safe run tag
function tag(; L, T, p)
    f(x) = replace(@sprintf("%.3f", x), "." => "p")
    return "L$(L)_T$(T)_p$(f(p))"
end

function main(args)
    opts = parse_args(args, build_parser())
    L    = opts["L"]
    T    = opts["T"]
    p    = opts["p"]

    observables = [:χ2]

    # Hardcoded output dir for easy local testing
    outdir = joinpath(ROOT, "output")
    mkpath(outdir)

    _, _ = circuit(neel_state(DiagonalStateMPS, 2), 2, 2, 0.1; observables=observables)

    t1 = time()
    state = neel_state(DiagonalStateMPS, L)
    state, data = circuit(state, L, T, p; observables=observables)
    dt = time() - t1

    tagstr  = tag(L = L, T = T, p = p)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))

    fname = joinpath(outdir, "sample_$(tagstr)_$(timestr)_$(randtag).jld2")

    @save fname L T p observables data dt
end

main(ARGS)
