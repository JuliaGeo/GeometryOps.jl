@inline Base.:-(x::CrossingFloat) = raw(-x.hi, -x.lo, x.rad)

@inline function add(x::CrossingFloat, y::CrossingFloat, check_domain::Val)
    sum_hi, sum_error = twosum(x.hi, y.hi)
    if (x.lo == 0.0) & (y.lo == 0.0) & (x.rad == 0.0) & (y.rad == 0.0)
        return finish(check_domain, sum_hi, sum_error, 0.0)
    end
    partial_low = sum_error + x.lo
    combined_low = partial_low + y.lo
    result_hi, result_lo = twosum(sum_hi, combined_low)
    error_bound = x.rad + y.rad + UNIT_ROUNDOFF * (
        (abs(sum_error) + abs(x.lo)) + (abs(partial_low) + abs(y.lo)))
    return finish(check_domain, result_hi, result_lo, pad(error_bound))
end

@inline Base.:+(x::CrossingFloat, y::CrossingFloat) = add(x, y, Val(true))
@inline Base.:-(x::CrossingFloat, y::CrossingFloat) = x + (-y)

@inline function mul(x::CrossingFloat, y::CrossingFloat, check_domain::Val)
    if x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0
        return isbounded(y) ? raw(0.0, 0.0, 0.0) : bad()
    elseif y.hi == 0.0 && y.lo == 0.0 && y.rad == 0.0
        return isbounded(x) ? raw(0.0, 0.0, 0.0) : bad()
    end
    product_hi, product_lo = twoprod(x.hi, y.hi)
    if (x.lo == 0.0) & (y.lo == 0.0) & (x.rad == 0.0) & (y.rad == 0.0)
        return finish(check_domain, product_hi, product_lo, 0.0)
    end
    high_low = x.hi * y.lo
    low_high = x.lo * y.hi
    cross_terms = high_low + low_high
    combined_low = product_lo + cross_terms
    result_hi, result_lo = twosum(product_hi, combined_low)
    rounding_error = UNIT_ROUNDOFF * (
        2 * (abs(high_low) + abs(low_high)) + abs(product_lo) + abs(cross_terms)) +
        abs(x.lo * y.lo)
    error_bound = rounding_error +
        x.rad * (abs(y.hi) + abs(y.lo) + y.rad) +
        (abs(x.hi) + abs(x.lo)) * y.rad
    return finish(check_domain, result_hi, result_lo, pad(error_bound))
end

@inline Base.:*(x::CrossingFloat, y::CrossingFloat) = mul(x, y, Val(true))

# Shared power-of-two scaling preserves x/y and puts its denominator in this domain.
@inline function _div_normalized(x::CrossingFloat, y::CrossingFloat)
    denominator_hi = y.hi
    denominator_lo = y.lo
    denominator_radius = y.rad
    ((1.0 <= abs(denominator_hi) <= 4.0) &&
     abs(denominator_lo) + denominator_radius <= 0.19 * abs(denominator_hi)) || return bad()
    x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0 && return raw(0.0, 0.0, 0.0)
    quotient_hi = x.hi / denominator_hi
    admissible(quotient_hi) || return bad()
    product_hi, product_lo = twoprod(quotient_hi, denominator_hi)
    residual_hi, residual_lo = twosum(x.hi, -product_hi)
    residual_after_product = residual_hi - product_lo
    residual_after_sum = residual_after_product + residual_lo
    residual_with_low = residual_after_sum + x.lo
    low_product = quotient_hi * denominator_lo
    correction_residual = residual_with_low - low_product
    quotient_correction = correction_residual / denominator_hi
    result_hi, result_lo = twosum(quotient_hi, quotient_correction)
    residual_rounding_error = UNIT_ROUNDOFF * (
        (abs(residual_hi) + abs(product_lo)) +
        (abs(residual_after_product) + abs(residual_lo)) +
        (abs(residual_after_sum) + abs(x.lo)) +
        abs(low_product) + (abs(residual_with_low) + abs(low_product)))
    # Center uncertainty and input radius contribute separate denominator bounds.
    center_lower = (abs(denominator_hi) - abs(denominator_lo)) * BOUND_DEFLATION
    interval_lower = (center_lower - denominator_radius) * BOUND_DEFLATION
    error_bound = UNIT_ROUNDOFF * abs(correction_residual) / abs(denominator_hi) +
        abs(correction_residual) * abs(denominator_lo) / (abs(denominator_hi) * center_lower) +
        residual_rounding_error / center_lower +
        x.rad / interval_lower +
        (abs(x.hi) + abs(x.lo) + x.rad) * denominator_radius / (center_lower * interval_lower)
    exact = (residual_hi == 0.0) & (residual_lo == 0.0) & (product_lo == 0.0) & (x.lo == 0.0) &
            (denominator_lo == 0.0) & (x.rad == 0.0) & (denominator_radius == 0.0)
    return make(result_hi, result_lo, exact ? 0.0 : pad(error_bound))
end

@inline function Base.:/(x::CrossingFloat, y::CrossingFloat)
    denominator_hi = y.hi
    1.0 <= abs(denominator_hi) <= 4.0 && return _div_normalized(x, y)
    isbounded(x) & isbounded(y) || return bad()
    denominator_hi == 0.0 && return bad()
    abs(y.lo) + y.rad <= 0.19 * abs(denominator_hi) || return bad()
    normalization = PowerOfTwo(-exponent(denominator_hi))
    scaled_numerator = scale_pow2(x, normalization)
    scaled_denominator = scale_pow2(y, normalization)
    isbounded(scaled_numerator) & isbounded(scaled_denominator) || return bad()
    return _div_normalized(scaled_numerator, scaled_denominator)
end

@inline function Base.sqrt(x::CrossingFloat)
    input_hi, input_lo, input_radius = x.hi, x.lo, x.rad
    input_hi == 0.0 && input_lo == 0.0 && input_radius == 0.0 && return raw(0.0, 0.0, 0.0)
    ((1.0 <= input_hi <= 16.0) && abs(input_lo) + input_radius <= 0.19 * input_hi) || return bad()
    root_hi = sqrt(input_hi)
    product_hi, product_lo = twoprod(root_hi, root_hi)
    residual_hi = input_hi - product_hi
    residual_after_product = residual_hi - product_lo
    correction_residual = residual_after_product + input_lo
    root_correction = correction_residual / (2 * root_hi)
    result_hi, result_lo = twosum(root_hi, root_correction)
    root_error_bound = abs(input_lo) / (1.8 * root_hi) + 2 * UNIT_ROUNDOFF * root_hi
    # One Newton correction leaves a quadratic error in the initial root estimate.
    newton_error = root_error_bound * root_error_bound / (2 * root_hi)
    rounding_error = (
        UNIT_ROUNDOFF * (abs(residual_hi) + abs(residual_after_product) + abs(input_lo)) +
        UNIT_ROUNDOFF * (abs(residual_hi) + abs(product_lo))) / (2 * root_hi) +
        UNIT_ROUNDOFF * abs(correction_residual) / (2 * root_hi)
    error_bound = newton_error + rounding_error + input_radius / (0.9 * root_hi)
    exact = (residual_hi == 0.0) & (product_lo == 0.0) & (input_lo == 0.0) & (input_radius == 0.0)
    return make(result_hi, result_lo, exact ? 0.0 : pad(error_bound))
end
