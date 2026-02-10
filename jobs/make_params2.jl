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
Ts      = 1000:1000:10000
gamma   = 0.1
Bs      = 0:5:100
samples = 10000

# -------------------------------
# Generate params.txt
# -------------------------------
lines = String[]
for T in Ts, B in Bs
    push!(lines, format_line(L=L, T=T, gamma=gamma, B=B, samples=samples))
end

write_params(lines; append=append)
