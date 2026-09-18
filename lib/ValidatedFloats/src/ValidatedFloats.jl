module ValidatedFloats

export ValidatedFloat, center, radius, isbounded, certify, certified_sign, CrossingFloats

include("CrossingFloats.jl")

"""
    ValidatedFloat(x::Float64)

A double-double approximation together with a conservative absolute error
radius. The represented real value lies within `radius(x)` of the exact sum of
the stored `hi` and `lo` limbs while
`isbounded(x)` is true. Arithmetic deliberately becomes unbounded when its
finite, normal-range preconditions cannot be established.
"""
struct ValidatedFloat <: Real
    hi::Float64
    lo::Float64
    rad::Float64
    function ValidatedFloat(hi::Float64, lo::Float64, rad::Float64)
        isnan(rad) && throw(ArgumentError("radius must not be NaN"))
        rad < 0 && throw(ArgumentError("radius must be nonnegative"))
        uppergap = isfinite(nextfloat(hi)) ? abs(nextfloat(hi) - hi) : abs(hi - prevfloat(hi))
        if isfinite(rad) && (!isfinite(hi) || !isfinite(lo) ||
                (hi == 0 ? lo != 0 : abs(lo) > max(uppergap, abs(hi - prevfloat(hi)))))
            throw(ArgumentError("center must be a finite normalized two-limb expansion"))
        end
        new(hi, lo, rad)
    end
    @inline ValidatedFloat(::Val{:trusted}, hi::Float64, lo::Float64, rad::Float64) =
        new(hi, lo, rad)
end

@inline _trusted(hi::Float64, lo::Float64, rad::Float64) =
    ValidatedFloat(Val(:trusted), hi, lo, rad)
@inline ValidatedFloat(x::Float64) = isfinite(x) ? _trusted(x, 0.0, 0.0) :
    _trusted(x, 0.0, Inf)
"The nearest `Float64` to the two-limb center (the radius is not included)."
center(x::ValidatedFloat) = x.hi + x.lo
radius(x::ValidatedFloat) = x.rad
isbounded(x::ValidatedFloat) = isfinite(x.hi) & isfinite(x.lo) & isfinite(x.rad)

const U = 0x1p-53
const MIN_NORMAL = floatmin(Float64)

@inline function _twosum(a::Float64, b::Float64)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end

@inline function _twoproduct(a::Float64, b::Float64)
    p = a * b
    return p, fma(a, b, -p)
end

@inline _mag(x::ValidatedFloat) = abs(x.hi) + abs(x.lo) + x.rad
@inline _mag_up(x::ValidatedFloat) = _up(_up(abs(x.hi) + abs(x.lo)) + x.rad)
@inline function _lower_abs(x::ValidatedFloat)
    d = prevfloat(prevfloat(abs(x.hi) - abs(x.lo)) - x.rad)
    return max(0.0, d)
end
@inline _normal_or_zero(x::Float64) = x == 0.0 || abs(x) >= MIN_NORMAL
@inline _finite_normal(x::ValidatedFloat) = isbounded(x) && _normal_or_zero(x.hi) && _normal_or_zero(x.lo)
@inline _up(x::Float64) = x == 0.0 ? 0.0 : (isfinite(x) ? nextfloat(x) : Inf)
@inline function _add_up(a::Float64, b::Float64)
    z = a + b
    isfinite(z) || return Inf
    z == 0.0 ? ((a == 0.0 && b == 0.0) ? 0.0 : nextfloat(0.0)) : nextfloat(z)
end
@inline function _mul_up(a::Float64, b::Float64)
    z = a * b
    isfinite(z) || return Inf
    z == 0.0 ? ((a == 0.0 || b == 0.0) ? 0.0 : nextfloat(0.0)) : nextfloat(z)
end
@inline function _div_up(a::Float64, b::Float64)
    z = a / b
    isfinite(z) || return Inf
    z == 0.0 ? (a == 0.0 ? 0.0 : nextfloat(0.0)) : nextfloat(z)
end
@inline _sum_up(xs...) = foldl(_add_up, xs; init=0.0)
@inline _roundbound(scale::Float64, n::Int) = _mul_up(nextfloat((n * U) / (1 - n * U)), scale)
@inline _unbounded(h::Float64 = NaN) = _trusted(h, 0.0, Inf)

@inline Base.:-(x::ValidatedFloat) = _trusted(-x.hi, -x.lo, x.rad)

@inline function Base.:+(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi + y.hi)
    sh, se = _twosum(x.hi, y.hi)
    t = se + x.lo
    t2 = t + y.lo
    rh, rl = _twosum(sh, t2)
    vals = (sh, se, t, t2, rh, rl)
    (all(isfinite, vals) && all(_normal_or_zero, vals)) || return _unbounded(rh)
    scale = _sum_up(abs(se), abs(x.lo), abs(t), abs(y.lo), abs(t2))
    rad = _sum_up(x.rad, y.rad, _roundbound(scale, 4))
    return _trusted(rh, rl, rad)
end

@inline Base.:-(x::ValidatedFloat, y::ValidatedFloat) = x + (-y)

@inline function Base.:*(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi * y.hi)
    (x.hi == 0 || y.hi == 0 || exponent(x.hi) + exponent(y.hi) >= -900) ||
        return _unbounded(x.hi * y.hi)
    (x.hi == 0 || y.hi == 0 || abs(x.hi) >= MIN_NORMAL / abs(y.hi)) || return _unbounded(x.hi * y.hi)
    ph, pl = _twoproduct(x.hi, y.hi)
    (ph != 0 || x.hi == 0 || y.hi == 0) || return _unbounded(ph)
    c1 = x.hi * y.lo
    c2 = x.lo * y.hi
    ((c1 != 0 || x.hi == 0 || y.lo == 0) &&
     (c2 != 0 || x.lo == 0 || y.hi == 0)) || return _unbounded(ph)
    c = c1 + c2
    t = pl + c
    rh, rl = _twosum(ph, t)
    vals = (ph, pl, c1, c2, c, t, rh, rl)
    (all(isfinite, vals) && all(_normal_or_zero, vals)) || return _unbounded(rh)
    dropped = _mul_up(abs(x.lo), abs(y.lo))
    prop = _add_up(_mul_up(x.rad, _mag_up(y)),
                   _mul_up(_add_up(abs(x.hi), abs(x.lo)), y.rad))
    scale = _sum_up(abs(c1), abs(c2), abs(c), abs(pl), abs(t))
    rad = _sum_up(prop, dropped, _roundbound(scale, 5))
    return _trusted(rh, rl, rad)
end

@inline function Base.:/(x::ValidatedFloat, y::ValidatedFloat)
    (_finite_normal(x) && _finite_normal(y)) || return _unbounded(x.hi / y.hi)
    ymin = _lower_abs(y)
    (isfinite(ymin) && ymin > 0) || return _unbounded(x.hi / y.hi)
    q1 = x.hi / y.hi
    ph, pl = _twoproduct(q1, y.hi)
    sh, sl = _twosum(x.hi, -ph)
    t1 = sh - pl; t2 = t1 + sl; t3 = t2 + x.lo
    m = q1 * y.lo; r = t3 - m; q2 = r / y.hi
    rh, rl = _twosum(q1, q2)
    vals = (q1, ph, pl, sh, sl, t1, t2, t3, m, r, q2, rh, rl)
    (all(isfinite, vals) && all(_normal_or_zero, vals)) || return _unbounded(rh)
    qv = _trusted(rh, rl, 0.0)
    residual = x - qv * y
    isbounded(residual) || return _unbounded(rh)
    # |x/y-q| = |x-q*y|/|y|. Both numerator and denominator are
    # outward bounds produced by this same abstraction.
    denom = _lower_abs(y)
    denom > 0 || return _unbounded(rh)
    return _trusted(rh, rl, _div_up(_mag_up(residual), denom))
end

@inline function Base.sqrt(x::ValidatedFloat)
    _finite_normal(x) || return _unbounded(sqrt(x.hi))
    x.hi == 0 && x.lo == 0 && x.rad == 0 && return _trusted(0.0, 0.0, 0.0)
    lower = prevfloat(prevfloat(x.hi - abs(x.lo)) - x.rad)
    lower >= 0 || throw(DomainError(x, "interval contains negative values"))
    lower == 0 && !(x.hi == 0 && x.lo == 0 && x.rad == 0) && return _unbounded(0.0)
    x.hi > 0 || return _unbounded()
    s0 = sqrt(x.hi)
    ph, pl = _twoproduct(s0, s0)
    r1 = x.hi - ph; r2 = r1 - pl; r = r2 + x.lo
    q = r / (2s0)
    rh, rl = _twosum(s0, q)
    vals = (s0, ph, pl, r1, r2, r, q, rh, rl)
    (all(isfinite, vals) && all(_normal_or_zero, vals)) || return _unbounded(rh)
    qv = _trusted(rh, rl, 0.0)
    residual = x - qv * qv
    isbounded(residual) || return _unbounded(rh)
    # |sqrt(x)-q| = |x-q^2|/(sqrt(x)+q); q alone is a conservative
    # positive lower bound for the denominator.
    denom = _lower_abs(qv)
    denom > 0 || return _unbounded(rh)
    return _trusted(rh, rl, _div_up(_mag_up(residual), denom))
end

@inline function Base.ldexp(x::ValidatedFloat, n::Integer)
    isbounded(x) || return _unbounded(ldexp(x.hi, n))
    h, l, r = ldexp(x.hi, n), ldexp(x.lo, n), ldexp(x.rad, n)
    ((x.hi == 0 || h != 0) && (x.lo == 0 || l != 0) && (x.rad == 0 || r != 0) &&
     all(isfinite, (h, l, r)) && all(_normal_or_zero, (h, l, r))) || return _unbounded(h)
    return _trusted(h, l, r)
end

"""
    certify(Float64, x::ValidatedFloat) -> Union{Float64,Nothing}

Return the unique round-to-nearest `Float64` for every real represented by
`x`, or `nothing` when this cannot be proven. Midpoint ties are intentionally
not certified.
"""
@inline function certify(::Type{Float64}, x::ValidatedFloat)
    isbounded(x) || return nothing
    x.rad == 0 && x.lo == 0 && return x.hi
    f, rem = _twosum(x.hi, x.lo)
    isfinite(f) || return nothing
    rem == 0 && x.rad == 0 && return f
    f == 0 && return nothing
    gap = min(abs(nextfloat(f) - f), abs(f - prevfloat(f)))
    (_up(abs(rem) + x.rad) < 0.5gap) ? f : nothing
end

"""
    certified_sign(x::ValidatedFloat) -> Union{Int,Nothing}

Return `-1`, `0`, or `1` when the entire represented interval has that sign;
return `nothing` when it contains or may cross zero.
"""
@inline function certified_sign(x::ValidatedFloat)
    isbounded(x) || return nothing
    x.hi == 0 && x.lo == 0 && x.rad == 0 && return 0
    _lower_abs(x) > 0 || return nothing
    return signbit(x.hi) ? -1 : 1
end

Base.show(io::IO, x::ValidatedFloat) = print(io, "ValidatedFloat(", x.hi, ", ", x.lo, ", ±", x.rad, ")")

end
