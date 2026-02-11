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

function one_bra(site::Index)
    return sum([dag(ITensors.state(site, "$i")) for i in 0:dim(site)-1])
end

function lognorm(state::DiagonalStateMPS)
    return lognorm(state.mps)
end

function lognorm(ψ::MPS)
    sites = siteinds(ψ)
    L = length(sites)

    mats = Vector{ITensor}(undef, L)
    for j in 1:L
        s = sites[j]
        bra = one_bra(s)
        mats[j] = bra * ψ[j]
    end

    v₀ = mats[end]
    logscale = 0.0
    for j in 1:L-1
        v₁ = mats[end-j] * v₀
        n = norm(v₁)
        logscale += log(n)
        v₀ = v₁ / n
    end

    x = logscale + log(Complex(Array(v₀)[]))
    if abs(imag(x)) < 1e-10
        x = real(x)
    end

    return x
end

function normalize!(state::DiagonalStateMPS)
    logscale = lognorm(state)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)

    scale = exp(logscale / L)
    for i in 1:L
        state.mps[i] /= scale
    end
    return state
end

function number_left(ψ::MPS)
    sites = siteinds(ψ)
    L = length(sites)

    mats = Vector{ITensor}(undef, L-1)
    for j in 2:L
        s = sites[j]
        bra = one_bra(s)
        mats[j-1] = bra * ψ[j]
    end

    v₀ = mats[end]
    logscale = 0.0
    for j in 1:L-2
        v₁ = mats[end-j] * v₀
        n = norm(v₁)
        logscale += log(n)
        v₀ = v₁ / n
    end

    ns = ψ[1] * v₀ * exp(logscale)

    return Array(ns, inds(ns)...)
end

function simple_circuit(state::DiagonalStateMPS, L::Int, T::Int, p::Float64; cutoff=1e-8, maxdim=200)
    SWAPn1 = decoherence_layer(state, SWAP, p, 1:2:L-1)
    SWAPn2 = decoherence_layer(state, SWAP, p, 2:2:L-1)

    for t in 1:T
        state = apply(SWAPn1, state; cutoff=cutoff, maxdim=maxdim)
        state = apply(SWAPn2, state; cutoff=cutoff, maxdim=maxdim)
        normalize!(state)
        truncate!(state; cutoff=cutoff, maxdim=maxdim)
    end

    return state
end

function number_reduce_left(state::DiagonalStateMPS, A::Int)
    state = deepcopy(state)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)


    dist = ITensor(1.0)
    dist *= ψ[1]
    i1 = only(inds(ψ[1],"Site"))
    for (counter, site) in enumerate(1:A-1)
        i2 = only(inds(ψ[site+1],"Site"))
        i3 = siteind("Qudit", 0; dim=counter+2, conserve_number=true, addtags="Left")

        T = ITensor(dag(i1), dag(i2), i3)
        for a in 1:counter+1, b in 1:2
            c = a + b - 1
            T[i1 => a, i2 => b, i3 => c] = 1.0
        end

        dist *= ψ[site+1]*T
        i1 = i3
    end

    return DiagonalStateMPS(MPS([dist, ψ[A+1:end]...]))
end

function number_reduce_right(state::DiagonalStateMPS, A::Int)
    state = deepcopy(state)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites)


    dist = ITensor(1.0)
    dist *= ψ[end]
    i1 = only(inds(ψ[end],"Site"))
    for (counter, site) in enumerate(1:A-1)
        i2 = only(inds(ψ[end-site],"Site"))
        i3 = siteind("Qudit", -1; dim=counter+2, conserve_number=true, addtags="Right")

        T = ITensor(dag(i1), dag(i2), i3)
        for a in 1:counter+1, b in 1:2
            c = a + b - 1
            T[i1 => a, i2 => b, i3 => c] = 1.0
        end

        dist *= ψ[end-site]*T
        i1 = i3
    end

    return DiagonalStateMPS(MPS([ψ[1:end-A]..., dist]))
end

function expval(state::DiagonalStateMPS, M::AbstractMatrix, pos::Int; cutoff=1E-8, maxdim=200, refs=0, conserve_qns=false)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Mψ = apply(op(M, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...), ψ; cutoff=cutoff, maxdim=maxdim)
    
    val = exp(lognorm(Mψ) - lognorm(ψ))
    return val
end

function measure(state::DiagonalStateMPS, M::AbstractMatrix, λ::Float64, pos::Int, m::Bool; cutoff=1E-8, maxdim=200, refs=0)
    ψ = state.mps
    sites = siteinds(ψ)
    L = length(sites) - refs
    M_width = Int(log2(size(M)[1]))

    Π = (I + (-1)^m * λ*M)/(sqrt(2*(1+λ^2)))
    g = op(Π*Π, [sites[mod1(pos+i,L)] for i in 0:M_width-1]...)

    ψ = apply(g, ψ; cutoff=cutoff, maxdim=maxdim)
    
    state = DiagonalStateMPS(ψ)
    truncate!(state; cutoff=cutoff, maxdim=maxdim)
    normalize!(state)

    return state
end

function sample_decoder(L::Int, T::Int, gamma::Float64, B::Int, samples::Int)
    p = exp(-gamma)*sinh(gamma)
    state = neel_state(DiagonalStateMPS, L; conserve_number=true)
    state = simple_circuit(state, L, T, p)

    state = number_reduce_left(state, L÷2)
    state = number_reduce_right(state, L÷2-B)

    total = 0.0
    total2 = 0.0
    for _ in 1:samples
        sampled_state = deepcopy(state)
        sampled_state, _, _ = measure(sampled_state, PauliZ, 1.0, 2:2+B-1; cutoff=1e-8, maxdim=200)
        ns = number_left(sampled_state.mps)
        prob = max(abs.(ns)...)
        total += prob
        total2 += prob^2
    end
    return total/samples, total2/samples
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
        
    _, _ = sample_decoder(4, 2, 0.1, 1, 2)

    t1 = time()
    prob1, prob2 = sample_decoder(L, T, gamma, B, samples)

    dt = time() - t1

    tagstr  = tag(L = L, T = T, gamma = gamma, B = B, samples = samples)
    timestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    randtag = string(rand(UInt32))

    fname = joinpath(outdir, "optimal_decoder_2_$(tagstr)_$(timestr)_$(randtag).jld2")

    @save fname L T gamma B samples prob1 prob2 dt
end

main(ARGS)
