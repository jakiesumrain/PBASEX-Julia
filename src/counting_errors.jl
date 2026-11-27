module CountingErrors

using LinearAlgebra: norm

# Integral of each piecewise-linear interpolation function over the window.
function window_integrals(energy::Vector{T}, window) where T
    window isa Union{Tuple,AbstractVector} && length(window) == 2 &&
        all(x -> x isa Real, window) ||
        throw(ArgumentError("Each peak window needs two numeric energy limits"))
    low, high = T.(window)
    first(energy) <= low < high <= last(energy) ||
        throw(ArgumentError("Peak windows must lie within the sampled energy range and have positive width"))
    integrals = zeros(T, length(energy))
    for j in 1:length(energy)-1
        left, right = max(low, energy[j]), min(high, energy[j+1])
        left >= right && continue
        mean_fraction = ((left - energy[j]) + (right - energy[j])) / (2 * (energy[j+1] - energy[j]))
        integrals[j] += (right - left) * (1 - mean_fraction)
        integrals[j+1] += (right - left) * mean_fraction
    end
    return (low, high), integrals
end

function summarize_counting_errors(result, basis, singular_filter, count_variance, alpha, windows)
    T = eltype(result.E)
    if !isempty(windows) && !all(>(zero(T)), diff(result.E))
        throw(ArgumentError("Peak integration requires a strictly increasing energy grid"))
    end
    # Distinct folded positions use disjoint raw pixels in both geometries.
    # A column of this scaled adjoint has norm equal to the observable's SE.
    folded_se = sqrt.(vec(count_variance))
    detector_response(singular_coordinates) =
        folded_se .* (basis.U * (singular_filter .* singular_coordinates))
    energy_basis = basis.frk .* transpose(basis.x ./ (2 * alpha))
    radial_count, energy_count = size(energy_basis)
    isotropic_vectors = @view basis.V[1:radial_count, :]
    spectrum_se = similar(result.IE)
    # Bound scratch space on a full detector; never form its covariance.
    for first_index in 1:16:energy_count
        indices = first_index:min(first_index + 15, energy_count)
        response = detector_response(isotropic_vectors' * view(energy_basis, :, indices))
        spectrum_se[indices] .= norm.(eachcol(response))
    end

    peaks = map(windows) do limits
        window, integrals = window_integrals(result.E, limits)
        radial_integrals = energy_basis * integrals
        yields = result.coefficients' * radial_integrals
        order_count = length(yields)
        singular_coordinates = zeros(T, length(singular_filter), order_count)
        for order_index in eachindex(yields)
            rows = (order_index-1)*radial_count+1:order_index*radial_count
            singular_coordinates[:, order_index] .= view(basis.V, rows, :)' * radial_integrals
        end
        response = detector_response(singular_coordinates)
        isotropic = view(response, :, 1)
        yield_se = norm(isotropic)
        betas = iszero(yields[1]) ? fill(T(NaN), order_count-1) : yields[2:end] ./ yields[1]
        beta_se = fill(T(NaN), length(betas))
        if yields[1] > yield_se
            for j in eachindex(betas)
                # Shared numerator and denominator noise cancels before squaring.
                beta_se[j] = norm(view(response, :, j+1) .- betas[j] .* isotropic) / yields[1]
            end
        end
        (; window, yield=yields[1], yield_se, betas, beta_se)
    end
    return (; result, spectrum_se, peaks)
end

end # module CountingErrors
