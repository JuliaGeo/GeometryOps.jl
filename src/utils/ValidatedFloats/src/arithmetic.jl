@inline Base.:-(x::ValidatedFloat) = _trusted(-x.hi, -x.lo, x.rad)

@inline function Base.:+(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi + y.hi)
    sum_hi, sum_error = _twosum(x.hi, y.hi)
    partial_low = sum_error + x.lo
    combined_low = partial_low + y.lo
    result_hi, result_lo = _twosum(sum_hi, combined_low)
    intermediates = (sum_hi, sum_error, partial_low, combined_low, result_hi, result_lo)
    (all(isfinite, intermediates) && all(_normal_or_zero, intermediates)) || return _unbounded(result_hi)
    scale = _sum_up(abs(sum_error), abs(x.lo), abs(partial_low), abs(y.lo), abs(combined_low))
    error_radius = _sum_up(x.rad, y.rad, _roundbound(scale, 4))
    return _trusted(result_hi, result_lo, error_radius)
end

@inline Base.:-(x::ValidatedFloat, y::ValidatedFloat) = x + (-y)

@inline function Base.:*(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi * y.hi)
    (x.hi == 0 || y.hi == 0 || exponent(x.hi) + exponent(y.hi) >= -900) ||
        return _unbounded(x.hi * y.hi)
    (x.hi == 0 || y.hi == 0 || abs(x.hi) >= MIN_NORMAL / abs(y.hi)) || return _unbounded(x.hi * y.hi)
    product_hi, product_lo = _twoproduct(x.hi, y.hi)
    (product_hi != 0 || x.hi == 0 || y.hi == 0) || return _unbounded(product_hi)
    high_low = x.hi * y.lo
    low_high = x.lo * y.hi
    ((high_low != 0 || x.hi == 0 || y.lo == 0) &&
     (low_high != 0 || x.lo == 0 || y.hi == 0)) || return _unbounded(product_hi)
    cross_terms = high_low + low_high
    combined_low = product_lo + cross_terms
    result_hi, result_lo = _twosum(product_hi, combined_low)
    intermediates = (product_hi, product_lo, high_low, low_high, cross_terms, combined_low, result_hi, result_lo)
    (all(isfinite, intermediates) && all(_normal_or_zero, intermediates)) || return _unbounded(result_hi)
    low_product_bound = _mul_up(abs(x.lo), abs(y.lo))
    input_error = _add_up(_mul_up(x.rad, _mag_up(y)),
                   _mul_up(_add_up(abs(x.hi), abs(x.lo)), y.rad))
    scale = _sum_up(abs(high_low), abs(low_high), abs(cross_terms), abs(product_lo), abs(combined_low))
    error_radius = _sum_up(input_error, low_product_bound, _roundbound(scale, 5))
    return _trusted(result_hi, result_lo, error_radius)
end

@inline function Base.:/(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi / y.hi)
    denominator_lower = _lower_abs(y)
    (isfinite(denominator_lower) && denominator_lower > 0) || return _unbounded(x.hi / y.hi)
    quotient_hi = x.hi / y.hi
    product_hi, product_lo = _twoproduct(quotient_hi, y.hi)
    residual_hi, residual_lo = _twosum(x.hi, -product_hi)
    residual_after_product = residual_hi - product_lo
    residual_after_sum = residual_after_product + residual_lo
    residual_with_low = residual_after_sum + x.lo
    low_product = quotient_hi * y.lo
    correction_residual = residual_with_low - low_product
    quotient_correction = correction_residual / y.hi
    result_hi, result_lo = _twosum(quotient_hi, quotient_correction)
    intermediates = (
        quotient_hi, product_hi, product_lo, residual_hi, residual_lo,
        residual_after_product, residual_after_sum, residual_with_low,
        low_product, correction_residual, quotient_correction, result_hi, result_lo)
    (all(isfinite, intermediates) && all(_normal_or_zero, intermediates)) || return _unbounded(result_hi)
    candidate = _trusted(result_hi, result_lo, 0.0)
    residual = x - candidate * y
    isbounded(residual) || return _unbounded(result_hi)
    # The quotient error is bounded by |x - candidate*y| / |y|.
    denominator_bound = _lower_abs(y)
    denominator_bound > 0 || return _unbounded(result_hi)
    return _trusted(result_hi, result_lo, _div_up(_mag_up(residual), denominator_bound))
end

@inline function Base.sqrt(x::ValidatedFloat)
    _finite_normal(x) || return _unbounded(sqrt(x.hi))
    x.hi == 0 && x.lo == 0 && x.rad == 0 && return _trusted(0.0, 0.0, 0.0)
    lower = prevfloat(prevfloat(x.hi - abs(x.lo)) - x.rad)
    lower >= 0 || throw(DomainError(x, "interval contains negative values"))
    lower == 0 && !(x.hi == 0 && x.lo == 0 && x.rad == 0) && return _unbounded(0.0)
    x.hi > 0 || return _unbounded()
    root_hi = sqrt(x.hi)
    product_hi, product_lo = _twoproduct(root_hi, root_hi)
    residual_hi = x.hi - product_hi
    residual_after_product = residual_hi - product_lo
    correction_residual = residual_after_product + x.lo
    root_correction = correction_residual / (2 * root_hi)
    result_hi, result_lo = _twosum(root_hi, root_correction)
    intermediates = (
        root_hi, product_hi, product_lo, residual_hi, residual_after_product,
        correction_residual, root_correction, result_hi, result_lo)
    (all(isfinite, intermediates) && all(_normal_or_zero, intermediates)) || return _unbounded(result_hi)
    candidate = _trusted(result_hi, result_lo, 0.0)
    residual = x - candidate * candidate
    isbounded(residual) || return _unbounded(result_hi)
    # The positive candidate bounds sqrt(x) + candidate from below.
    denominator_bound = _lower_abs(candidate)
    denominator_bound > 0 || return _unbounded(result_hi)
    return _trusted(result_hi, result_lo, _div_up(_mag_up(residual), denominator_bound))
end

@inline function Base.ldexp(x::ValidatedFloat, n::Integer)
    isbounded(x) || return _unbounded(ldexp(x.hi, n))
    scaled_hi, scaled_lo, scaled_radius = ldexp(x.hi, n), ldexp(x.lo, n), ldexp(x.rad, n)
    ((x.hi == 0 || scaled_hi != 0) && (x.lo == 0 || scaled_lo != 0) && (x.rad == 0 || scaled_radius != 0) &&
     all(isfinite, (scaled_hi, scaled_lo, scaled_radius)) && all(_normal_or_zero, (scaled_hi, scaled_lo, scaled_radius))) || return _unbounded(scaled_hi)
    return _trusted(scaled_hi, scaled_lo, scaled_radius)
end
