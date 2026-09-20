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
        upper_gap = isfinite(nextfloat(hi)) ? abs(nextfloat(hi) - hi) : 0.0
        lower_gap = isfinite(prevfloat(hi)) ? abs(hi - prevfloat(hi)) : 0.0
        if isfinite(rad) && (!isfinite(hi) || !isfinite(lo) ||
                (hi == 0 ? lo != 0 : abs(lo) > max(upper_gap, lower_gap)))
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

"""
    certify(Float64, x::ValidatedFloat) -> Union{Float64, Nothing}

Return the unique round-to-nearest `Float64` for every real represented by
`x`, or `nothing` when this cannot be proven. Midpoint ties are intentionally
not certified.
"""
@inline function certify(::Type{Float64}, x::ValidatedFloat)
    isbounded(x) || return nothing
    x.rad == 0 && x.lo == 0 && return x.hi
    rounded, remainder = _twosum(x.hi, x.lo)
    isfinite(rounded) || return nothing
    remainder == 0 && x.rad == 0 && return rounded
    rounded == 0 && return nothing
    # The smaller neighbor gap gives a safe rounding cell on either side.
    neighbor_gap = min(abs(nextfloat(rounded) - rounded), abs(rounded - prevfloat(rounded)))
    (_up(abs(remainder) + x.rad) < 0.5 * neighbor_gap) ? rounded : nothing
end

"""
    certified_sign(x::ValidatedFloat) -> Union{Int, Nothing}

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
