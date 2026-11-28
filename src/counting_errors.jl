module CountingErrors

using LinearAlgebra: mul!, norm

# Integral of each piecewise-linear interpolation function over the window.
function window_integrals(energy::Vector{T}, window) where T
    integrals = zeros(T, length(energy))
    return window_integrals!(integrals, energy, window), integrals
end

function window_integrals!(integrals, energy::Vector{T}, window) where T
    window isa Union{Tuple,AbstractVector} && length(window) == 2 &&
        all(x -> x isa Real, window) ||
        throw(ArgumentError("Each peak window needs two numeric energy limits"))
    low, high = T.(window)
    first(energy) <= low < high <= last(energy) ||
        throw(ArgumentError("Peak windows must lie within the sampled energy range and have positive width"))
    fill!(integrals, zero(T))
    for j in 1:length(energy)-1
        left, right = max(low, energy[j]), min(high, energy[j+1])
        left >= right && continue
        mean_fraction = ((left - energy[j]) + (right - energy[j])) / (2 * (energy[j+1] - energy[j]))
        integrals[j] += (right - left) * (1 - mean_fraction)
        integrals[j+1] += (right - left) * mean_fraction
    end
    return low, high
end

function summarize_counting_errors(result, basis, singular_filter, count_variance, alpha, windows;
                                   native_errors=true)
    T = eltype(result.E)
    if !isempty(windows) && !all(>(zero(T)), diff(result.E))
        throw(ArgumentError("Peak integration requires a strictly increasing energy grid"))
    end
    # Distinct folded positions use disjoint raw pixels in both geometries.
    # A column of this scaled adjoint has norm equal to the observable's SE.
    folded_se = sqrt.(vec(count_variance))
    energy_basis = basis.frk .* transpose(basis.x ./ (2 * alpha))
    spectrum_se, beta_se = if native_errors
        native_counting_errors(result, basis, singular_filter, folded_se, energy_basis)
    else
        nothing, nothing
    end
    peaks = window_counting_errors(result, basis, singular_filter, folded_se, energy_basis, windows)
    return (; result, spectrum_se, peaks, beta_se, windows=peaks, legendre_orders=result.legendre_orders)
end

function native_counting_errors(result, basis, singular_filter, folded_se, energy_basis)
    T = eltype(result.E)
    detector_response(singular_coordinates) =
        folded_se .* (basis.U * (singular_filter .* singular_coordinates))
    radial_count, energy_count = size(energy_basis)
    isotropic_vectors = @view basis.V[1:radial_count, :]
    spectrum_se = similar(result.IE)
    beta_se = fill(T(NaN), size(result.betas))
    # Bound scratch space on a full detector; never form its covariance.
    for first_index in 1:16:energy_count
        indices = first_index:min(first_index + 15, energy_count)
        radial_samples = view(energy_basis, :, indices)
        response = detector_response(isotropic_vectors' * radial_samples)
        spectrum_se[indices] .= norm.(eachcol(response))
        for order in axes(beta_se, 2)
            rows = order*radial_count+1:(order+1)*radial_count
            angular_response = detector_response(view(basis.V, rows, :)' * radial_samples)
            for (column, i) in enumerate(indices)
                if result.IE[i] > spectrum_se[i]
                    beta_se[i, order] = ratio_se(view(angular_response, :, column),
                        view(response, :, column), result.betas[i, order], result.IE[i])
                end
            end
        end
    end
    return spectrum_se, beta_se
end

function window_counting_errors(result, basis, singular_filter, folded_se, energy_basis, windows)
    T = eltype(result.E)
    peak_type = NamedTuple{(:window, :yield, :yield_se, :betas, :beta_se, :spectrum, :spectrum_se),
                           Tuple{Tuple{T,T},T,T,Vector{T},Vector{T},T,T}}
    peaks = similar(windows, peak_type)
    isempty(windows) && return peaks
    radial_count, energy_count = size(energy_basis)
    order_count = size(result.coefficients, 2)
    batch_size = min(8, length(windows))
    limits = Vector{Tuple{T,T}}(undef, batch_size)
    integrals = Vector{T}(undef, energy_count)
    radial = Matrix{T}(undef, radial_count, batch_size)
    yields = Matrix{T}(undef, order_count, batch_size)
    coordinates = Matrix{T}(undef, length(singular_filter), order_count * batch_size)
    responses = Matrix{T}(undef, length(folded_se), order_count * batch_size)
    for first_window in 1:batch_size:length(windows)
        count = min(batch_size, length(windows) - first_window + 1)
        for column in 1:count
            limits[column] = window_integrals!(integrals, result.E, windows[first_window + column - 1])
            radial_column = view(radial, :, column)
            mul!(radial_column, energy_basis, integrals)
            mul!(view(yields, :, column), result.coefficients', radial_column)
        end
        radial_batch = view(radial, :, 1:count)
        for order in 1:order_count
            rows = (order-1)*radial_count+1:order*radial_count
            cols = (order-1)*count+1:order*count
            mul!(view(coordinates, :, cols), view(basis.V, rows, :)', radial_batch)
        end
        active = 1:order_count*count
        coordinates_batch = view(coordinates, :, active)
        coordinates_batch .*= singular_filter
        response = view(responses, :, active)
        mul!(response, basis.U, coordinates_batch)
        response .*= folded_se
        for column in 1:count
            window = limits[column]
            yield = yields[1, column]
            isotropic = view(response, :, column)
            yield_se = norm(isotropic)
            betas = iszero(yield) ? fill(T(NaN), order_count-1) : yields[2:end, column] ./ yield
            beta_se = fill(T(NaN), length(betas))
            if yield > yield_se
                for j in eachindex(betas)
                    angular = view(response, :, j*count + column)
                    # Keep numerator and denominator noise from the same pixels.
                    angular .-= betas[j] .* isotropic
                    beta_se[j] = norm(angular) / yield
                end
            end
            width = window[2] - window[1]
            peaks[first_window + column - 1] =
                (; window, yield, yield_se, betas, beta_se, spectrum=yield/width, spectrum_se=yield_se/width)
        end
    end
    return peaks
end

# Subtract the shared numerator/denominator noise before taking the norm.
ratio_se(angular, isotropic, beta, intensity) = norm(angular .- beta .* isotropic) / intensity

end # module CountingErrors
