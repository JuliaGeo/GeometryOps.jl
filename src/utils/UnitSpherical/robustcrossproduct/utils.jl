
# Helper function to find a perpendicular vector to a given vector
# This is a generic implementation that works with any vector type
function find_orthogonal(v::AbstractVector)
    # Choose the smallest component to zero out
    x, y, z = v[1], v[2], v[3]
    ax, ay, az = abs(x), abs(y), abs(z)
    
    if ax <= ay && ax <= az
        # x is smallest, zero it out
        return UnitSphericalPoint(0.0, -z, y)
    elseif ay <= ax && ay <= az
        # y is smallest, zero it out
        return UnitSphericalPoint(-z, 0.0, x)
    else
        # z is smallest, zero it out
        return UnitSphericalPoint(-y, x, 0.0)
    end
end

"""
    isUnitLength(v::AbstractVector)

Check if a vector has unit length within a small tolerance.

Returns true if the vector's magnitude is approximately 1.0.
"""
function isUnitLength(v::AbstractVector)
    return isapprox(sum(abs2, v), 1.0, rtol=1e-14)
end

"""
    isNormalizable(v::AbstractVector)

Returns true if the given vector's magnitude is large enough such that
the angle to another vector of the same magnitude can be measured using 
angle calculations without loss of precision due to floating-point underflow.

This matches S2's IsNormalizable function.
"""
function isNormalizable(v::AbstractVector)
    # Same threshold as in S2 - the largest component should be at least 2^-242
    # This ensures we can normalize without precision loss
    return maximum(abs.(v)) >= ldexp(1.0, -242)
end

"""
    normalization_needed(v::AbstractVector)

Determines if a vector's magnitude is too small for reliable normalization.
Returns true if the vector needs special handling to avoid numerical issues.

This is essentially the opposite of isNormalizable.
"""
function normalization_needed(v::AbstractVector)
    # The threshold is 1e-14 squared, approximately 1e-28
    norm_v² = sum(abs2, v)
    return norm_v² < 1e-28
end

"""
    ensureNormalizable(p::AbstractVector)

Scales a 3-vector as necessary to ensure that the result can be normalized
without loss of precision due to floating-point underflow.

This matches S2's EnsureNormalizable function.
"""
function ensureNormalizable(p::AbstractVector)
    if p == zeros(eltype(p), 3)
        error("Vector must be non-zero")
    end
    
    if !isNormalizable(p)
        # Scale so that the largest component has magnitude in [1, 2)
        p_max = maximum(abs.(p))
        # Scale by 2^(-1-exponent) to achieve this range
        return ldexp(2.0, -1 - exponent(Float64(p_max))) * p
    end
    
    return p
end

"""
    normalizableFromExact(xf::Vector{BigFloat})

Converts a BigFloat vector to a double-precision vector, scaling the
result as necessary to ensure that the result can be normalized without 
loss of precision due to floating-point underflow.

This matches S2's NormalizableFromExact function.
"""
function normalizableFromExact(xf::AbstractVector{BigFloat})
    # First try a simple conversion
    x = Float64.(xf)
    
    if isNormalizable(x)
        return x
    end
    
    # Find the largest exponent
    max_exp = -1000000  # Very small initial value
    for i in 1:3
        if !iszero(xf[i])
            max_exp = max(max_exp, exponent(xf[i]))
        end
    end
    
    if max_exp < -1000000  # No non-zero components
        return zero(xf)
    end
    
    # Scale to get components in a good range
    return Float64.(ldexp.(Float64.(xf), -max_exp))
end


"""
    isless_vector(a::AbstractVector, b::AbstractVector)

Lexicographic comparison of vectors.
This is used to establish a consistent order for symbolic perturbations.

Returns true if a comes before b in lexicographic order.
"""
function isless_vector(a::AbstractVector, b::AbstractVector)
    # First compare x coordinates
    a[1] != b[1] && return a[1] < b[1]
    
    # If x coordinates are equal, compare y coordinates
    a[2] != b[2] && return a[2] < b[2]
    
    # If both x and y are equal, compare z coordinates
    return a[3] < b[3]
end


# Use the Base.< function for UnitSphericalPoint to delegate to our isless_vector function
# TODO: this is sailing the high seas!
function Base.:<(a::UnitSphericalPoint, b::UnitSphericalPoint)
    return isless_vector(a, b)
end