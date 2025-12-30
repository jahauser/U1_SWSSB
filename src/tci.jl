using ITensorMPS, ITensors, LinearAlgebra, Random#,  QuanticsTCI
import TensorCrossInterpolation as TCI
import TensorCrossInterpolation: GlobalPivotSearchInput, MultiIndex, crossinterpolate2

struct SlaterGlobalPivotFinder <: TCI.AbstractGlobalPivotFinder
    nsearch::Int
    maxnglobalpivot::Int
    tolmarginglobalsearch::Float64
    sysNp::Int
end

function SlaterGlobalPivotFinder(sysNP::Int;
    nsearch::Int=5,
    maxnglobalpivot::Int=5,
    tolmarginglobalsearch::Float64=10.0
)
    return SlaterGlobalPivotFinder(nsearch, maxnglobalpivot, tolmarginglobalsearch, sysNP)
end

function (finder::SlaterGlobalPivotFinder)(
    input::GlobalPivotSearchInput{ValueType},
    f,
    abstol::Float64;
    verbosity::Int=0,
    rng::AbstractRNG=Random.default_rng()
)::Vector{MultiIndex} where {ValueType}
    L = length(input.localdims)
    nsearch = finder.nsearch
    maxnglobalpivot = finder.maxnglobalpivot
    tolmarginglobalsearch = finder.tolmarginglobalsearch

    # Generate random initial points (with defined particle number)
    initial_points = [random_binary_vector(L, finder.sysNp; inds=false) .+ 1 for _ in 1:nsearch]

    # Find pivots with high interpolation error
    found_pivots = MultiIndex[]
    for point in initial_points
        # Perform local search from each initial point
        current_point = copy(point)
        best_error = 0.0
        best_point = copy(point)

        # Local search
        for p in 1:L
            for v in 1:input.localdims[p]
                current_point[p] = v
                error = abs(f(current_point) - input.current_tt(current_point))
                if error > best_error
                    best_error = error
                    best_point = copy(current_point)
                end
            end
            current_point[p] = point[p]  # Reset to original point
        end

        # Add point if error is above threshold
        if best_error > abstol * tolmarginglobalsearch
            push!(found_pivots, best_point)
        end
    end

    # Limit number of pivots
    if length(found_pivots) > maxnglobalpivot
        found_pivots = found_pivots[1:maxnglobalpivot]
    end

    if verbosity > 0
        println("Found $(length(found_pivots)) global pivots")
    end

    return found_pivots
end

function get_prob(ψ0::MPS, config)
  N = length(ψ0)
  sites = siteinds(ψ0)
  amplitude = ITensor(1.0)
  for j in 1:N
    sⱼ = sites[j]
    Aⱼ = ψ0[j]
    stateⱼ = ITensorMPS.state(sⱼ, config[j])
    amplitude *= Aⱼ * dag(stateⱼ)
  end
  return abs((amplitude)[])
end

function get_prob_vec(ψ0::MPS, v)
  config = [(iseven(vi) ? "1" : "0") for vi in v]
  return get_prob(ψ0, config)
end

function get_sqrtprob_vec(ψ0::MPS, v)
  config = [(iseven(vi) ? "1" : "0") for vi in v]
  return sqrt(get_prob(ψ0, config))
end

function random_binary_vector(n, m; inds=true)
  v = zeros(Int, n)
  indices = randperm(n)[1:m]
  if inds
    return indices
  end
  v[indices] .= 1
  return v
end

function tci(state::DiagonalStateMPS; n_initialpivots=40, tol=1e-8, maxdim=3*maxlinkdim(state.mps))
    L = length(state.mps)
    # state = normalize(state)
    ψ = state.mps
    normalize!(ψ)
    orthogonalize!(ψ, 1)

    f(v) = get_sqrtprob_vec(ψ, v)
    localdims = fill(2, L)
    initialpivots = [sample(ψ) for _ in 1:n_initialpivots]
    tci, ranks, errors = TCI.crossinterpolate2(Float64, f, localdims, initialpivots; tolerance=tol, maxbonddim=maxdim,
                                               globalpivotfinder=SlaterGlobalPivotFinder(L÷2))


    s = siteinds("Qubit", L; conserve_qns=false)
    ψ1 = MPS(tci; sites=s)
    return ψ1/exp(lognorm(ψ1))
end