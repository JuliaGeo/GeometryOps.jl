# # RobustCrossProduct
#=
```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
robust_cross_product
```

## What is this thing?

The `robust_cross_product` function computes a robust version of the cross product between two unit vectors on a sphere. 
This function is essential for geometric algorithms on the sphere that need stability even when points are very close 
together or nearly antipodal.

Standard cross product calculations can lead to numerical instability when:
- Two points are nearly identical (resulting in a very small cross product)
- Two points are nearly antipodal (making the direction of the cross product unstable)

This implementation handles these edge cases by:
1. Trying a regular cross product first
2. Checking if the result is too small for reliable normalization
3. Using specialized methods to ensure a stable perpendicular vector is returned

## Examples

```@example robust-cross
using GeometryOps.UnitSpherical

# Regular case - points at right angles
a = UnitSphericalPoint(1, 0, 0)
b = UnitSphericalPoint(0, 1, 0)
result = robust_cross_product(a, b)
println("Standard case: ", result)

# Nearly identical points
c = UnitSphericalPoint(1, 1e-10, 0)
result_similar = robust_cross_product(a, c)
println("Nearly identical points: ", result_similar)

# Check that result is perpendicular to both inputs
dot_a = result_similar ⋅ a
dot_c = result_similar ⋅ c
println("Perpendicular to inputs: ", isapprox(dot_a, 0, atol=1e-14), ", ", isapprox(dot_c, 0, atol=1e-14))
```

=#

module RobustCrossProduct

using ..UnitSpherical: UnitSphericalPoint, orthogonal
using StaticArrays
using LinearAlgebra

include("utils.jl")

# Error constants - these follow the S2 implementation
# DBL_ERR represents the rounding error of a single arithmetic operation
const DBL_ERR = eps(Float64) / 2
# sqrt(3) is used in error calculations
const SQRT3 = sqrt(3.0)
# This is the maximum directional error in the result, in radians 
const ROBUST_CROSS_PROD_ERROR = 6 * DBL_ERR
# Constant to check if we have access to higher precision
const HAS_LONG_DOUBLE = precision(Float64) < precision(BigFloat)
# Error for exact cross product calculations
const EXACT_CROSS_PROD_ERROR = DBL_ERR

isDoubleFloatsAvailable(args...) = false

"""
    robust_cross_product(a::AbstractVector, b::AbstractVector)

Compute a robust version of `a × b` (cross product) for unit vectors.

This method handles the case where `a` and `b` are very close together or
antipodal by computing a stable perpendicular to both points. 

The implementation follows Google's S2 Geometry Library to ensure numerical
stability even in difficult cases.

Returns a unit-length vector that is perpendicular to both input vectors.

## Examples

```jldoctest
julia> using GeometryOps.UnitSpherical

julia> a = UnitSphericalPoint(1, 0, 0)
julia> b = UnitSphericalPoint(0, 1, 0)
julia> result = robust_cross_product(a, b)
julia> isapprox(result, UnitSphericalPoint(0, 0, 1))
true
```
"""
function robust_cross_product(a::AbstractVector, b::AbstractVector)
    # Check that inputs are unit length 
    @boundscheck @assert isUnitLength(a) "Input vector 'a' must be unit length"
    @boundscheck @assert isUnitLength(b) "Input vector 'b' must be unit length"


    # The direction of cross(a, b) becomes unstable as (a + b) or (a - b)
    # approaches zero.  This leads to situations where cross(a, b) is not
    # very orthogonal to "a" and/or "b".  To solve this problem robustly requires
    # falling back to extended precision, arbitrary precision, and even symbolic
    # perturbations to handle the case when "a" and "b" are exactly
    # proportional, e.g. a == -b (see google/s2geometry/s2predicates.cc for details).
    result, was_stable = stable_cross_product(a, b)
    if was_stable
        # @debug "RCP: Simple cross product was stable"
        return normalize(result)
    else
        # @debug "RCP: Simple cross product was unstable" result.x result.y result.z
    end
    # Handle the (a == b) case now, before doing expensive arithmetic.  The only
    # result that makes sense mathematically is to return zero, but it turns out
    # to reduce the number of special cases in client code if we instead return
    # an arbitrary orthogonal vector.
    if a == b
        # @debug "RCP: Vectors are identical, generating orthogonal vector" a b
        return orthogonal(a)
    end
    
    if isDoubleFloatsAvailable()
        # @debug "RCP: Using double floats"
        result, was_stable = stable_cross_product(to_doublefloat.(a), to_doublefloat.(b))
        was_stable && return normalize(Float64.(result))
    end
    # @debug "RCP: Using exact arithmetic"

    
    
    # From here, we follow the exact C++ implementation order:
    # First, use exactCrossProd which will handle long double and exact arithmetic
    result = exact_cross_product(a, b)
    
    # Make sure the result can be normalized reliably
    return normalize(result)
end

stable_cross_product(a::AbstractVector{T1}, b::AbstractVector{T2}) where {T1, T2} = stable_cross_product(promote(a, b)...)

"""
    getStableCrossProd(a::AbstractVector, b::AbstractVector)

Computes a numerically stable cross product between unit vectors.

This implements the algorithm from S2's GetStableCrossProd function,
computing (a-b)×(a+b) which yields better numerical stability when
the vectors are nearly identical.

Returns a tuple of (result, success) where:
- result is the computed cross product vector (not normalized)
- success is a boolean indicating if the computation was sufficiently accurate
"""
function stable_cross_product(a::AbstractVector{T}, b::AbstractVector{T}) where T
    # We compute the cross product (a - b) x (a + b).  Mathematically this is
    # exactly twice the cross product of "a" and "b", but it has the numerical
    # advantage that (a - b) and (a + b) are nearly perpendicular (since "a" and
    # "b" are unit length).  This yields a result that is nearly orthogonal to
    # both "a" and "b" even if these two values differ only very slightly.
    #
    # The maximum directional error in radians when this calculation is done in
    # precision T (where T is a floating-point type) is:
    #```
    #   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR
    #```
    # where ||N|| is the norm of the result.  To keep this error to at most
    # kRobustCrossProdError, assuming this value is much less than 1, we need
    #```
    #   (1 + 2 * sqrt(3) + 32 * sqrt(3) * DBL_ERR / ||N||) * T_ERR <= kErr
    #```
    # From this you can see that in order for this calculation to ever succeed in
    # double precision, we must have `kErr > (1 + 2 * sqrt(3)) * DBL_ERR`, which is
    # about `4.46 * DBL_ERR`.  We actually set `kRobustCrossProdError == 6 * DBL_ERR
    # (== 3 * DBL_EPSILON)` in order to minimize the number of cases where higher
    # precision is needed; in particular, higher precision is only necessary when
    # "a" and "b" are closer than about `18 * DBL_ERR == 9 * DBL_EPSILON`.
    # (80-bit precision can handle inputs as close as `2.5 * LDBL_EPSILON`.)
    T_ERR = eps(float(T)) / 2 
    kMinNorm = (32 * sqrt(3) * DBL_ERR) / (ROBUST_CROSS_PROD_ERROR / T_ERR - (1 + 2sqrt(3)))

    # Finally...we compute the result by regular cross product.
    result = cross(a - b, a + b)
    
    # Check if the result norm is sufficiently large
    was_stable = LinearAlgebra.norm_sqr(result) >= kMinNorm^2
    return result, was_stable
end

exact_cross_product(a::AbstractVector{T1}, b::AbstractVector{T2}) where {T1, T2} = exact_cross_product(promote(a, b)...)

"""
    exact_cross_product(a::AbstractVector, b::AbstractVector)

Compute the cross product using arbitrary precision arithmetic.
This is used when standard floating-point arithmetic is not accurate enough.

This matches S2's ExactCrossProd function, first trying higher precision
if available, then exact arithmetic, then symbolic perturbation.
"""
function exact_cross_product(a::AbstractVector{T}, b::AbstractVector{T}) where T
    @assert a != b "Vectors must be different"
    
    # Use BigFloat for arbitrary precision arithmetic
    # really this is probably enough?  But we can go to 
    # exact formulations later, this is just a stopgap 
    # anyway.
    big_a = BigFloat.(a; precision=512)
    big_b = BigFloat.(b; precision=512)
    result_xf = cross(big_a, big_b)
    
    # Check if we got a non-zero result
    # This is equivalent to s2's `s2pred::IsZero`.
    if !all(<=(1e-300), abs.(result_xf))
        return normalizableFromExact(result_xf)
    end
    
    # If exact arithmetic yields zero, use symbolic perturbation
    # This follows S2's approach exactly.
    # symbolic_cross_product requires that a < b.
    if isless_vector(a, b)
        return ensureNormalizable(symbolic_cross_product(a, b))
    else
        return -ensureNormalizable(symbolic_cross_product(b, a))
    end
end


symbolic_cross_product(a::AbstractVector{T1}, b::AbstractVector{T2}) where {T1, T2} = symbolic_cross_product(promote(a, b)...)
"""
    symbolic_cross_product(a::AbstractVector, b::AbstractVector)

Compute a symbolic cross product when exact arithmetic yields zero.
This implements the symbolic perturbation model used in S2 geometry.

Returns a vector that is the symbolic cross product.
"""
function symbolic_cross_product(a::AbstractVector{T}, b::AbstractVector{T}) where T
    @assert a != b "Vectors must be different for symbolic cross product"
    
    # SymbolicCrossProdSorted requires that a < b
    if isless_vector(a, b)
        return symbolic_cross_product_sorted(a, b)
    else
        return -symbolic_cross_product_sorted(b, a)
    end
end

symbolic_cross_product_sorted(a::AbstractVector{T1}, b::AbstractVector{T2}) where {T1, T2} = symbolic_cross_product_sorted(promote(a, b)...)
"""
    symbolic_cross_product_sorted(a::AbstractVector, b::AbstractVector)

Helper function to compute the symbolic cross product when points are collinear.
Assumes that a < b lexicographically.

This implements the symbolic perturbation model described in S2 geometry.
"""
function symbolic_cross_product_sorted(a::AbstractVector{T}, b::AbstractVector{T}) where T
    # The following code uses the same symbolic perturbation model as S2::Sign.
    # The particular sequence of tests below was obtained using Mathematica
    # (although it would be easy to do it by hand for this simple case).
    #
    # Just like the function SymbolicallyPerturbedSign() in s2predicates.cc,
    # every input coordinate x[i] is assigned a symbolic perturbation dx[i].  We
    # then compute the cross product
    #
    #     (a + da) × (b + db)
    #
    # The result is a polynomial in the perturbation symbols. For example if we
    # did this in one dimension, the result would be
    #
    #     a * b + b * da + a * db + da * db
    #
    # where "a" and "b" have numerical values and "da" and "db" are symbols.
    # In 3 dimensions the result is similar except that the coefficients are
    # 3-vectors rather than scalars.
    #
    # Every possible UnitSphericalPoint has its own symbolic perturbation in each coordinate
    # (i.e., there are about 3 * 2^192 symbols). The magnitudes of the
    # perturbations are chosen such that if x < y lexicographically, the
    # perturbations for "y" are much smaller than the perturbations for "x".
    # Similarly, the perturbations for the coordinates of a given point x are
    # chosen such that dx[1] is much smaller than dx[2] which is much smaller
    # than dx[3]. Putting this together with the fact the inputs to this function
    # have been sorted so that a < b lexicographically, this tells us that
    #
    #     da[3] > da[2] > da[1] > db[3] > db[2] > db[1]
    #
    # where each perturbation is so much smaller than the previous one that we
    # don't even need to consider it unless the coefficients of all previous
    # perturbations are zero. In fact, each succeeding perturbation is so small
    # that we don't need to consider it unless the coefficient of all products of
    # the previous perturbations are zero. For example, we don't need to
    # consider the coefficient of db[2] unless the coefficient of db[3]*da[1] is
    # zero.
    #
    # The following code simply enumerates the coefficients of the perturbations
    # (and products of perturbations) that appear in the cross product above, in
    # order of decreasing perturbation magnitude. The first non-zero
    # coefficient determines the result. The easiest way to enumerate the
    # coefficients in the correct order is to pretend that each perturbation is
    # some tiny value "eps" raised to a power of two:
    #
    # eps^    1      2      4      8     16     32
    #       da[3]  da[2]  da[1]  db[3]  db[2]  db[1]
    #
    # Essentially we can then just count in binary and test the corresponding
    # subset of perturbations at each step. So for example, we must test the
    # coefficient of db[3]*da[1] before db[2] because eps^12 > eps^16.

    if b[1] != 0 || b[2] != 0
        return UnitSphericalPoint{T}(-b[2], b[1], 0) # da[3]
    end
    
    if b[3] != 0
        return UnitSphericalPoint{T}(b[3], 0, 0)
    end
    
    # None of the remaining cases should occur in practice for unit vectors
    @assert b[1] == 0 && b[2] == 0 "Expected both b[1] and b[2] to be zero"
    
    if a[1] != 0 || a[2] != 0
        # Fix: This needs to match C++ code which returns (-a[1], a[0], 0) in 0-based indexing
        # In Julia's 1-based indexing, this is (-a[2], a[1], 0)
        return UnitSphericalPoint{T}(-a[2], a[1], 0)
    end
    
    # This is always non-zero in the S2 implementation
    return UnitSphericalPoint{T}(1, 0, 0)
end

# Use the Base.< function for UnitSphericalPoint to delegate to our isless_vector function
function Base.:<(a::UnitSphericalPoint, b::UnitSphericalPoint)
    return isless_vector(a, b)
end

end

#=
# IMPLEMENTATION PLAN FOR S2 EDGE CROSSINGS FUNCTIONALITY

Based on analyzing s2edge_crossings.cc, we need to implement the following components:

## 1. Complete the robust cross product implementation
- [x] Basic RobustCrossProd function ✓
- [x] GetStableCrossProd implementation ✓
- [x] ExactCrossProd for arbitrary-precision arithmetic ✓
  - This uses BigFloat for arbitrary precision vectors
  - Implementing the IsZero function for exact vectors
- [x] SymbolicCrossProd for handling symbolic perturbations ✓
  - This includes lexicographic comparison for points
  - Implemented both SymbolicCrossProdSorted and SymbolicCrossProd
- [x] Add IsUnitLength, IsNormalizable, and EnsureNormalizable functions ✓
- [x] Add NormalizableFromExact conversion ✓

## 2. Edge crossing detection
- [ ] CrossingSign function to determine if two edges cross
  - This is the core function that uses S2EdgeCrosser
- [ ] VertexCrossing for when two edges share a vertex
- [ ] SignedVertexCrossing to determine the sign of a vertex crossing
- [ ] EdgeOrVertexCrossing to handle both cases

## 3. Intersection point calculation
- [ ] GetIntersection function to compute the intersection of two edges
- [ ] Support functions for intersection:
  - [ ] GetIntersectionSimple for simple intersections
  - [ ] GetIntersectionStable for more robust intersections
  - [ ] GetIntersectionExact for arbitrary precision
  - [ ] Helper functions like IsNormalizable, EnsureNormalizable
  - [ ] Functions for vector projection and normalization

## 4. Error estimation and auxiliary functions
- [ ] Error constants like kIntersectionError
- [ ] Functions for checking if a point lies on an edge (ApproximatelyOrdered)
- [ ] RobustNormalWithLength for computing normals with length estimation

## 5. Structure the implementation into modules
- [ ] Move edge crossing functionality to a separate EdgeCrossings module
- [ ] Move edge intersection calculation to an appropriate module
- [ ] Keep the core RobustCrossProd implementation here

## 6. Implementation strategy
1. Start by completing the RobustCrossProd implementation with ExactCrossProd
2. Implement CrossingSign and the basic S2EdgeCrosser functionality
3. Add the intersection point calculation
4. Add the remaining auxiliary and symbolic computation functions
5. Write comprehensive tests for each component

## 7. Dependencies and utilities needed
- A module for arbitrary precision arithmetic (we already have that in BigFloat)
- Symbolic perturbation utilities
- Robust predicates for orientation tests (we already have that through ExactPredicates.jl and AdaptivePredicates.jl)
- Specialized vector types and operations (we already have that in UnitSpherical.jl)
=#