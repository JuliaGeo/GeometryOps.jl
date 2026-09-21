# A restricted arithmetic domain keeps the leading-product residual representable.
const MIN_MAGNITUDE = 0x1p-400
const MAX_MAGNITUDE = 0x1p400
const UNIT_ROUNDOFF = 0x1p-53
const BOUND_INFLATION = 1.0 + 0x1p-44
const BOUND_DEFLATION = 1.0 - 0x1p-44
const UNDERFLOW_ALLOWANCE = 0x1p-1000

struct CrossingFloat <: Real
    hi::Float64
    lo::Float64
    rad::Float64
    CrossingFloat(hi::Float64, lo::Float64, error_radius::Float64, ::Val{:raw}) = new(hi, lo, error_radius)
end

@inline raw(hi, lo, error_radius) = CrossingFloat(hi, lo, error_radius, Val(:raw))
@inline bad() = raw(0.0, 0.0, Inf)
@inline finish(::Val{true}, hi, lo, error_radius) = make(hi, lo, error_radius)
@inline finish(::Val{false}, hi, lo, error_radius) = raw(hi, lo, error_radius)
@inline admit(x::CrossingFloat) = make(x.hi, x.lo, x.rad)
@inline function admissible(hi)
    # Positive Float64 bit patterns are ordered; unsigned subtraction rejects values below the range.
    bits = reinterpret(UInt64, hi) & 0x7fffffffffffffff
    lower_bits = reinterpret(UInt64, MIN_MAGNITUDE)
    upper_bits = reinterpret(UInt64, MAX_MAGNITUDE)
    return (bits == 0) | ((bits - lower_bits) <= (upper_bits - lower_bits))
end

@inline function make(hi, lo, error_radius)
    return admissible(hi) & (error_radius <= MAX_MAGNITUDE) ? raw(hi, lo, error_radius) : bad()
end

@inline CrossingFloat(x::Float64) = make(x, 0.0, 0.0)
@inline center(x::CrossingFloat) = x.hi + x.lo
@inline radius(x::CrossingFloat) = x.rad
@inline isbounded(x::CrossingFloat) = isfinite(x.rad)

# Inflation covers bound arithmetic; the additive allowance covers underflow.
@inline pad(error_radius) = (error_radius + UNDERFLOW_ALLOWANCE) * BOUND_INFLATION
@inline function twosum(a, b)
    sum_hi = a + b
    recovered_b = sum_hi - a
    return sum_hi, (a - (sum_hi - recovered_b)) + (b - recovered_b)
end

@inline function twoprod(a, b)
    product_hi = a * b
    return product_hi, fma(a, b, -product_hi)
end

struct PowerOfTwo
    value::Float64
    function PowerOfTwo(n::Integer)
        scale = ldexp(1.0, n)
        isfinite(scale) && scale > 0.0 || throw(ArgumentError("unrepresentable power of two"))
        new(scale)
    end
end

@inline function scale_pow2(x::CrossingFloat, power::PowerOfTwo)
    scale = power.value
    scaled_hi, scaled_lo, scaled_radius = x.hi * scale, x.lo * scale, x.rad * scale
    # Normal scaled limbs preserve the exact power-of-two transformation.
    ((x.lo == 0.0 || abs(scaled_lo) >= floatmin(Float64)) &&
     (x.rad == 0.0 || scaled_radius >= floatmin(Float64)) && (x.hi == 0.0 || scaled_hi != 0.0)) || return bad()
    return make(scaled_hi, scaled_lo, scaled_radius)
end

@inline Base.ldexp(x::CrossingFloat, n::Integer) = scale_pow2(x, PowerOfTwo(n))

@inline function certified_sign(x::CrossingFloat)
    isbounded(x) || return nothing
    x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0 && return 0
    abs(x.hi) > (abs(x.lo) + x.rad) * BOUND_INFLATION || return nothing
    return signbit(x.hi) ? -1 : 1
end

@inline function certify(::Type{Float64}, x::CrossingFloat)
    isbounded(x) || return nothing
    rounded, remainder = twosum(x.hi, x.lo)
    remainder == 0.0 && x.rad == 0.0 && return rounded
    rounded == 0.0 && return nothing
    # The smaller neighbor gap gives a safe rounding cell on either side.
    neighbor_gap = abs(rounded) - prevfloat(abs(rounded))
    (abs(remainder) + x.rad) * BOUND_INFLATION < 0.5 * neighbor_gap || return nothing
    return rounded
end
