const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
function format_line(; L, T, p, mode, t=nothing, samples=nothing)
    parts = String[
        "--L $L",
        "--T $T",
        "--p $p",
        "--mode $mode",
    ]
    if t !== nothing
        push!(parts, "--t $t")
    end
    if samples !== nothing
        push!(parts, "--samples $samples")
    end
    return join(parts, " ")
end

"Write or append parameter lines to jobs/params.txt."
function write_params(lines::Vector{String}; append::Bool=false)
    open(PARAMS_FILE, append ? "a" : "w") do io
        for ln in lines
            println(io, ln)
        end
    end
    action = append ? "Appended" : "Wrote"
    println("$action $(length(lines)) lines to $(PARAMS_FILE)")
end

# -------------------------------
# Sweeps
# -------------------------------
p = 0.1
append = false

T_of_L(L) = 10L  # ALWAYS

# L ranges
Ls_renyi2        = 4 .* (2 .^ (0:8))  # 4 .. 1024
Ls_renyi1_others = 4 .* (2 .^ (0:5))  # 4 .. 128

ts_renyi2(L) = [1]


ts_renyi1_TCI(L) = 0:L:T_of_L(L)
ts_renyi1_sampled(L) = 0:L:T_of_L(L)

# samples list acts as repeats, as a FUNCTION of L
# e.g. fill(100, 10) => 10 repeats with --samples 100
function samples_renyi1_sampled(L)
    if L ≤ 32
        return fill(5000, 1)
    elseif L == 64
        return fill(1000, 5)
    elseif L == 128
        return fill(100, 10)
    else
        return Int[]
    end
end

# -------------------------------
# Generate params.txt
# -------------------------------
lines = String[]

# renyi2: 9 jobs
for L in Ls_renyi2
    for t in ts_renyi2(L)  # typically [1]
        push!(lines, format_line(L=L, T=T_of_L(L), p=p, mode="renyi2", t=t))
    end
end

# renyi1_TCI: 60 jobs
for L in Ls_renyi1_others
    for t in ts_renyi1_TCI(L)
        push!(lines, format_line(L=L, T=T_of_L(L), p=p, mode="renyi1_TCI", t=t))
    end
end

# renyi1_sampled (samples list acts as repeats): 60*16=960
# for L in Ls_renyi1_others
#     ts = ts_renyi1_sampled(L)
#     ss = samples_renyi1_sampled(L)
#     for t in ts, s in ss
#         push!(lines, format_line(L=L, T=T_of_L(L), p=p, mode="renyi1_sampled", t=t, samples=s))
#     end
# end

write_params(lines; append=append)
