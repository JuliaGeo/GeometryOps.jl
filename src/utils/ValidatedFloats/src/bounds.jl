const UNIT_ROUNDOFF = 0x1p-53
const MIN_NORMAL = floatmin(Float64)

@inline function _twosum(a::Float64, b::Float64)
    sum_hi = a + b
    recovered_b = sum_hi - a
    return sum_hi, (a - (sum_hi - recovered_b)) + (b - recovered_b)
end

@inline function _twoproduct(a::Float64, b::Float64)
    product_hi = a * b
    return product_hi, fma(a, b, -product_hi)
end

@inline _mag(x::ValidatedFloat) = abs(x.hi) + abs(x.lo) + x.rad
@inline _mag_up(x::ValidatedFloat) = _up(_up(abs(x.hi) + abs(x.lo)) + x.rad)
@inline function _lower_abs(x::ValidatedFloat)
    lower_bound = prevfloat(prevfloat(abs(x.hi) - abs(x.lo)) - x.rad)
    return max(0.0, lower_bound)
end

@inline _normal_or_zero(x::Float64) = x == 0.0 || abs(x) >= MIN_NORMAL
@inline _finite_normal(x::ValidatedFloat) = isbounded(x) && _normal_or_zero(x.hi) && _normal_or_zero(x.lo)
@inline _up(x::Float64) = x == 0.0 ? 0.0 : (isfinite(x) ? nextfloat(x) : Inf)
# Positive underflow still needs a positive upper bound.
@inline function _add_up(a::Float64, b::Float64)
    rounded = a + b
    isfinite(rounded) || return Inf
    rounded == 0.0 ? ((a == 0.0 && b == 0.0) ? 0.0 : nextfloat(0.0)) : nextfloat(rounded)
end

@inline function _mul_up(a::Float64, b::Float64)
    rounded = a * b
    isfinite(rounded) || return Inf
    rounded == 0.0 ? ((a == 0.0 || b == 0.0) ? 0.0 : nextfloat(0.0)) : nextfloat(rounded)
end

@inline function _div_up(a::Float64, b::Float64)
    rounded = a / b
    isfinite(rounded) || return Inf
    rounded == 0.0 ? (a == 0.0 ? 0.0 : nextfloat(0.0)) : nextfloat(rounded)
end

@inline _sum_up(xs...) = foldl(_add_up, xs; init=0.0)
@inline _roundbound(scale::Float64, n::Int) = _mul_up(nextfloat((n * UNIT_ROUNDOFF) / (1 - n * UNIT_ROUNDOFF)), scale)
@inline _unbounded(h::Float64 = NaN) = _trusted(h, 0.0, Inf)
