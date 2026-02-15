const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
function format_line(; L, T, gamma, B, samples=nothing)
    parts = String[
        "--L $L",
        "--T $T",
        "--gamma $gamma",
        "--B $B",
    ]
    if samples !== nothing
        push!(parts, "--samples $samples")
    end
    return join(parts, " ")
end

"Write or append parameter lines to params.txt."
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
# Sweep
# -------------------------------
append  = false

L       = 200
Ts      = union(10:10:100)
gamma   = 0.1
Bs      = union(0:5:95, [99])
samples = 1000

num_repeats = Dict(
    0 => 8,
    5 => 8,
    10 => 8,
    15 => 8,
    20 => 8,
    25 => 8,
    30 => 8,
    35 => 8,
    40 => 8,
    45 => 8,
    50 => 8,
    55 => 8,
    60 => 8,
    65 => 8,
    70 => 8,
    75 => 8,
    80 => 8,
    85 => 8,
    90 => 8,
    95 => 8,
    99 => 10, # more repeats for harder cases
)

# -------------------------------
# Generate params.txt
# -------------------------------
lines = String[]

for T in Ts, B in Bs
    for repeats in 1:num_repeats[B]
        push!(lines, format_line(L=L, T=T, gamma=gamma, B=B, samples=samples))
    end
end

write_params(lines; append=append)
