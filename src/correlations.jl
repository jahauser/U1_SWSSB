# function correlation_sweep(ψ::MPS, op1::String, op2::String)
#     ψ /= exp(lognorm(ψ))
#     L = length(ψ)

#     output_sweep = zeros(ComplexF64, L÷2)

#     middle_stacks = Vector{ITensor}(undef, L÷2)
#     middle_stacks[1] = ITensor(1)
#     for k in 1:L÷2-1
#         i = L÷2 - k + 1
#         j = L÷2 + k
#         si = siteinds(ψ)[i]
#         sj = siteinds(ψ)[j]
#         middle_stacks[k+1] = middle_stacks[k] * ψ[i] * dag(prime(ψ[i], !si)) * ψ[j] * dag(prime(ψ[j], !sj))
#     end

#     left_stack = ITensor(1)
#     right_stack = ITensor(1)
#     for i in 1:L÷2
#         j = L - i + 1

#         si = siteinds(ψ)[i]
#         sj = siteinds(ψ)[j]

#         mat_i = ψ[i] * op(op1, siteinds(ψ)[i]) * dag(prime(ψ[i]))
#         mat_j = ψ[j] * op(op2, siteinds(ψ)[j]) * dag(prime(ψ[j]))
        
#         output_sweep[L÷2-i+1] = reduce(*, [left_stack, mat_i, middle_stacks[L÷2-i+1], mat_j, right_stack])[]
#         left_stack = left_stack * ψ[i] * dag(prime(ψ[i], !si))
#         right_stack = ψ[j] * dag(prime(ψ[j], !sj)) * right_stack
#     end

#     return output_sweep
# end