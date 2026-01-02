#=
# Minimum Bounding Circle

The minimum bounding circle (also called smallest enclosing circle) of a set of points
is the smallest circle that contains all the points.

GeometryOps provides Welzl's algorithm for computing this in expected O(n) time.

## Example

```julia
import GeometryOps as GO, GeoInterface as GI

points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
circle = GO.minimum_bounding_circle(points)
GO.radius(circle)  # approximately 0.707
```
=#

export minimum_bounding_circle, PlanarCircle, radius

"""
    PlanarCircle{T}

A circle in 2D Euclidean space, represented by its center and squared radius.

Use `radius(circle)` to get the actual radius (computes sqrt).

!!! warning "Experimental"
    This type is not part of the public API and is subject to change without a breaking version.
    It implements GeoInterface's `PolygonTrait`, so code using GeoInterface methods will remain
    compatible even if the concrete type changes.

## Fields
- `center::Tuple{T, T}`: The (x, y) coordinates of the center
- `radius_squared::T`: The squared radius (avoids sqrt in distance comparisons)
"""
struct PlanarCircle{T}
    center::Tuple{T, T}
    radius_squared::T
end

"""
    radius(circle::PlanarCircle)

Return the radius of the circle. Computes `sqrt(radius_squared)`.
"""
radius(c::PlanarCircle) = sqrt(c.radius_squared)

"""
    PlanarCircleRing{T}

A lazy wrapper that presents a `PlanarCircle` as a `LinearRing` with interpolated points.
Used internally for GeoInterface compatibility.

!!! warning "Internal"
    This is an internal type and not part of the public API.
"""
struct PlanarCircleRing{T}
    circle::PlanarCircle{T}
end

# GeoInterface for PlanarCircle (PolygonTrait)
GI.geomtrait(::PlanarCircle) = GI.PolygonTrait()
GI.nring(::GI.PolygonTrait, ::PlanarCircle) = 1
GI.getring(::GI.PolygonTrait, c::PlanarCircle, i::Int) = PlanarCircleRing(c)
GI.getexterior(::GI.PolygonTrait, c::PlanarCircle) = PlanarCircleRing(c)
GI.gethole(::GI.PolygonTrait, ::PlanarCircle) = ()

function GI.npoint(::GI.PolygonTrait, c::PlanarCircle)
    return GI.npoint(GI.LinearRingTrait(), PlanarCircleRing(c))
end

# GeoInterface for PlanarCircleRing (LinearRingTrait)
GI.geomtrait(::PlanarCircleRing) = GI.LinearRingTrait()
GI.npoint(::GI.LinearRingTrait, ::PlanarCircleRing) = 101  # 100 segments + closing point

function GI.getpoint(::GI.LinearRingTrait, r::PlanarCircleRing, i::Int)
    n = GI.npoint(GI.LinearRingTrait(), r)
    nsegs = n - 1
    idx = mod1(i, nsegs)
    θ = 2π * (idx - 1) / nsegs
    rad = radius(r.circle)
    sinθ, cosθ = sincos(θ) # faster than sin/cos separately
    x = r.circle.center[1] + rad * cosθ
    y = r.circle.center[2] + rad * sinθ
    return (x, y)
end

"""
    Welzl{M <: Manifold} <: ManifoldIndependentAlgorithm{M}

Welzl's algorithm for computing the minimum bounding circle.

This is a randomized algorithm with expected O(n) time complexity.
Works on any manifold given an appropriate distance function.

## Constructor

    Welzl(; manifold=Planar())

## Example

```julia
import GeometryOps as GO

points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
circle = GO.minimum_bounding_circle(GO.Welzl(), points)
```
"""
struct Welzl{M <: Manifold} <: ManifoldIndependentAlgorithm{M}
    manifold::M
end

Welzl(; manifold::Manifold=Planar()) = Welzl(manifold)

GeometryOpsCore.manifold(alg::Welzl) = alg.manifold

# Helper: check if point is inside or on circle (using squared distance)
function _point_in_circle(::Planar, p, circle::PlanarCircle)
    isnan(circle.radius_squared) && return false
    dx = p[1] - circle.center[1]
    dy = p[2] - circle.center[2]
    dist_squared = dx * dx + dy * dy
    # Use small epsilon for floating point tolerance
    return dist_squared <= circle.radius_squared * (1 + eps(Float64) * 10)
end

# Create circle with diameter from p1 to p2
function _circle_from_two_points(m::Planar, p1, p2)
    T = promote_type(typeof(p1[1]), typeof(p2[1]))
    center = ((p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)
    dx = p1[1] - center[1]
    dy = p1[2] - center[2]
    radius_squared = dx * dx + dy * dy
    return PlanarCircle{T}(center, radius_squared)
end

# Create circumcircle of three points
function _circle_from_three_points(m::Planar, p1, p2, p3)
    T = promote_type(typeof(p1[1]), typeof(p2[1]), typeof(p3[1]))

    ax, ay = p1[1], p1[2]
    bx, by = p2[1], p2[2]
    cx, cy = p3[1], p3[2]

    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))

    # Collinear points - fall back to two-point circle using farthest pair
    if abs(d) < eps(T) * 100
        d12 = (ax - bx)^2 + (ay - by)^2
        d23 = (bx - cx)^2 + (by - cy)^2
        d13 = (ax - cx)^2 + (ay - cy)^2
        if d12 >= d23 && d12 >= d13
            return _circle_from_two_points(m, p1, p2)
        elseif d23 >= d13
            return _circle_from_two_points(m, p2, p3)
        else
            return _circle_from_two_points(m, p1, p3)
        end
    end

    ux = ((ax^2 + ay^2) * (by - cy) + (bx^2 + by^2) * (cy - ay) + (cx^2 + cy^2) * (ay - by)) / d
    uy = ((ax^2 + ay^2) * (cx - bx) + (bx^2 + by^2) * (ax - cx) + (cx^2 + cy^2) * (bx - ax)) / d

    center = (ux, uy)
    radius_squared = (ax - ux)^2 + (ay - uy)^2

    return PlanarCircle{T}(center, radius_squared)
end

# Create minimum circle from 0-3 boundary points
function _make_circle(m::Planar, boundary::Vector)
    if isempty(boundary)
        return PlanarCircle((NaN, NaN), NaN)
    elseif length(boundary) == 1
        T = typeof(boundary[1][1])
        return PlanarCircle{T}(boundary[1], zero(T))
    elseif length(boundary) == 2
        return _circle_from_two_points(m, boundary[1], boundary[2])
    else
        return _circle_from_three_points(m, boundary[1], boundary[2], boundary[3])
    end
end

# Recursive Welzl algorithm
# points: all points to consider
# idx: current index (1-based, processes points[idx:end])
# boundary: points that must be on the circle boundary (0-3 points)
function _welzl!(m::Manifold, points::Vector, idx::Int, boundary::Vector)
    # Base case: no more points or 3 boundary points define a unique circle
    if idx > length(points) || length(boundary) == 3
        return _make_circle(m, boundary)
    end

    p = points[idx]

    # Recursively compute circle without p
    circle = _welzl!(m, points, idx + 1, boundary)

    # If p is inside the circle, we're done
    if _point_in_circle(m, p, circle)
        return circle
    end

    # Otherwise, p must be on the boundary of the minimum circle
    push!(boundary, p)
    result = _welzl!(m, points, idx + 1, boundary)
    pop!(boundary)

    return result
end

"""
    minimum_bounding_circle([algorithm], geometry)

Compute the minimum bounding circle of `geometry`.

Returns a circle geometry containing all points of the input. For planar geometries,
returns a `PlanarCircle`; for spherical geometries, returns a `SphericalCap`.

!!! warning "Return type subject to change"
    The concrete return type (currently `PlanarCircle` for planar manifold) may change
    in future versions without a breaking release. However, the return type will always
    implement GeoInterface, so code using GeoInterface methods (e.g., `GI.getexterior`,
    `GI.getpoint`) will remain compatible.

## Arguments
- `algorithm`: The algorithm to use. Defaults to `Welzl()` which uses Welzl's expected O(n) algorithm.
- `geometry`: Any geometry compatible with GeoInterface, or a vector of point-like objects.

## Example

```julia
import GeometryOps as GO, GeoInterface as GI

# From points
points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
circle = GO.minimum_bounding_circle(points)

# From any geometry
polygon = GI.Polygon([[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]])
circle = GO.minimum_bounding_circle(polygon)

# Access via GeoInterface for forward compatibility
ring = GI.getexterior(circle)
```
"""
function minimum_bounding_circle end

minimum_bounding_circle(geom) = minimum_bounding_circle(Welzl(), geom)

function minimum_bounding_circle(alg::Welzl{Planar}, geom)
    # Extract all points as tuples
    points = collect(flatten(tuples, GI.PointTrait, geom))

    # Handle edge cases
    if isempty(points)
        return PlanarCircle((NaN, NaN), NaN)
    elseif length(points) == 1
        T = typeof(points[1][1])
        return PlanarCircle{T}(points[1], zero(T))
    end

    # Shuffle for expected O(n) performance
    shuffled = Random.shuffle(points)

    # Initialize empty boundary
    T = typeof(shuffled[1][1])
    boundary = Tuple{T, T}[]

    return _welzl!(alg.manifold, shuffled, 1, boundary)
end


function minimum_bounding_circle(alg::Welzl{<: Spherical}, geom)
    # Extract all points as UnitSphericalPoints
    points = collect(flatten(UnitSpherical.UnitSphereFromGeographic(), GI.PointTrait, geom))

    # Handle edge cases
    if isempty(points)
        return SphericalCap(UnitSphericalPoint(NaN, NaN, NaN), NaN, NaN)
    elseif length(points) == 1
        T = typeof(points[1][1])
        return SphericalCap{T}(points[1], zero(T), one(T))
    end

    # Shuffle for expected O(n) performance
    shuffled = Random.shuffle(points)

    # Initialize empty boundary
    T = typeof(shuffled[1][1])
    boundary = UnitSphericalPoint{T}[]

    return _welzl!(alg.manifold, shuffled, 1, boundary)
end


# Create minimum circle from 0-3 boundary points
function _make_circle(m::Spherical, boundary::Vector)
    if isempty(boundary)
        return SphericalCap(UnitSphericalPoint(NaN, NaN, NaN), NaN, NaN)
    elseif length(boundary) == 1
        T = typeof(boundary[1][1])
        return SphericalCap{T}(boundary[1], zero(T), one(T))
    elseif length(boundary) == 2
        return _circle_from_two_points(m, boundary[1], boundary[2])
    else
        return _circle_from_three_points(m, boundary[1], boundary[2], boundary[3])
    end
end

function _circle_from_two_points(m::Spherical, p1::UnitSphericalPoint, p2::UnitSphericalPoint)
    midpoint = UnitSpherical.slerp(p1, p2, 0.5)
    radius = UnitSpherical.spherical_distance(p1, p2) / 2
    return SphericalCap(midpoint, radius)
end

function _circle_from_three_points(m::Spherical, p1::UnitSphericalPoint, p2::UnitSphericalPoint, p3::UnitSphericalPoint)
    return SphericalCap(p1, p2, p3)
end

# TODO: replace this internal calculation with some other thing using `distancelike`
function _point_in_circle(m::Spherical, p::UnitSphericalPoint, c::SphericalCap)
    isnan(c.radiuslike) && return false
    # For point inside cap: angular_distance(p, center) <= radius
    # Equivalently: cos(angular_distance) >= cos(radius)
    # Since p ⋅ center = cos(angular_distance) and radiuslike = cos(radius):
    return (p ⋅ c.point) >= c.radiuslike
end