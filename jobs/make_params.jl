const PARAMS_FILE = joinpath(@__DIR__, "params.txt")

"Format one CLI line for run.jl."
format_line(L, T, p) =
    "--L $L --T $T --p $p"

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
# Sweep
# -------------------------------
Ls = 16 .* (2 .^ (0:7))      # 16, 32, 64, ..., 2048
p = 0.1
append = false

lines = String[]
for L in Ls
    T = 2L
    push!(lines, format_line(L, T, p))
end

write_params(lines; append=append)
