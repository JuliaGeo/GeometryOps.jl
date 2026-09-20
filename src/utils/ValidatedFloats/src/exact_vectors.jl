# An admitted vector of exact Float64 coordinates. The tokenized inner
# constructor keeps unchecked construction inside this module.
struct _Exact3Token end
const _EXACT3_TOKEN = _Exact3Token()
struct Exact3
    xyz::NTuple{3, Float64}
    Exact3(xyz::NTuple{3, Float64}, ::_Exact3Token) = new(xyz)
end

"""
    exact3(xyz::NTuple{3, Float64}) -> Union{Exact3, Nothing}

Admit three exact Float64 coordinates once for the specialized vector kernels.
"""
@inline function exact3(xyz::NTuple{3, Float64})
    admissible(xyz[1]) & admissible(xyz[2]) & admissible(xyz[3]) || return nothing
    return Exact3(xyz, _EXACT3_TOKEN)
end

@inline function _exact_difference_product(a::Float64, b::Float64,
                                           c::Float64, d::Float64)
    product_hi, product_lo = twoprod(a, b)
    second_product_hi, second_product_lo = twoprod(c, d)
    (product_hi == second_product_hi) & (product_lo == second_product_lo) && return raw(0.0, 0.0, 0.0)
    difference_hi, difference_error = twosum(product_hi, -second_product_hi)
    partial_low = difference_error + product_lo
    combined_low = partial_low - second_product_lo
    result_hi, result_lo = twosum(difference_hi, combined_low)
    bound = UNIT_ROUNDOFF * ((abs(difference_error) + abs(product_lo)) + (abs(partial_low) + abs(second_product_lo)))
    return raw(result_hi, result_lo, pad(bound))
end

# Exact3 inputs make all input low limbs and radii statically zero. Each
# component pays one final admission after its two products are combined.
@inline function cross3(a::Exact3, b::Exact3)
    x = a.xyz
    y = b.xyz
    return (admit(_exact_difference_product(x[2], y[3], x[3], y[2])),
            admit(_exact_difference_product(x[3], y[1], x[1], y[3])),
            admit(_exact_difference_product(x[1], y[2], x[2], y[1])))
end

@inline function sum3(a::Exact3, b::Exact3)
    x = a.xyz
    y = b.xyz
    h1, l1 = twosum(x[1], y[1])
    h2, l2 = twosum(x[2], y[2])
    h3, l3 = twosum(x[3], y[3])
    return (make(h1, l1, 0.0), make(h2, l2, 0.0), make(h3, l3, 0.0))
end
