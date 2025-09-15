# # Spherical Caps

#=
```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
SphericalCap
circumcenter_on_unit_sphere
```

## What is SphericalCap?

A spherical cap represents a section of a unit sphere about some point, bounded by a radius. 
It is defined by a center point on the unit sphere and a radius (in radians).

Spherical caps are used in:
- Representing circular regions on a spherical surface
- Approximating and bounding spherical geometries
- Spatial indexing and filtering on the unit sphere
- Implementing containment, intersection, and disjoint predicates

The `SphericalCap` type offers multiple constructors to create caps from:
- UnitSphericalPoint and radius
- Geographic coordinates and radius
- Three points on the unit sphere (circumcircle)

## Examples

```@example sphericalcap
using GeometryOps
using GeoInterface

# Create a spherical cap from a point and radius
point = UnitSphericalPoint(1.0, 0.0, 0.0)  # Point on the unit sphere
cap = SphericalCap(point, 0.5)  # Cap with radius 0.5 radians
```

```@example sphericalcap
# Create a spherical cap from geographic coordinates
lat, lon = 40.0, -74.0  # New York City (approximate)
point = GeoInterface.Point(lon, lat)
cap = SphericalCap(point, 0.1)  # Cap with radius ~0.1 radians
```

```@example sphericalcap
# Create a spherical cap from three points (circumcircle)
p1 = UnitSphericalPoint(1.0, 0.0, 0.0)
p2 = UnitSphericalPoint(0.0, 1.0, 0.0)
p3 = UnitSphericalPoint(0.0, 0.0, 1.0)
cap = SphericalCap(p1, p2, p3)
```

=#

# Spherical cap implementation
struct SphericalCap{T}
    point::UnitSphericalPoint{T}
    radius::T
end

SphericalCap(point::UnitSphericalPoint{T}, radius::Number) where T = SphericalCap{T}(point, convert(T, radius))
SphericalCap(point, radius::Number) = SphericalCap(GI.trait(point), point, radius)
function SphericalCap(::GI.PointTrait, point, radius::Number)
    return SphericalCap(UnitSphereFromGeographic()(point), radius)
end

SphericalCap(geom) = SphericalCap(GI.trait(geom), geom)
SphericalCap(t::GI.PointTrait, geom) = SphericalCap(t, geom, 0)
# TODO: add implementations for line string and polygon traits
# TODO: add implementations for multitraits based on this

# TODO: this returns an approximately antipodal point...

# TODO: exact-predicate intersection
# This is all inexact and thus subject to floating point error
function _intersects(x::SphericalCap, y::SphericalCap)
    spherical_distance(x.point, y.point) <= x.radius + y.radius
end

_disjoint(x::SphericalCap, y::SphericalCap) = !_intersects(x, y)

function _contains(big::SphericalCap, small::SphericalCap)
    dist = spherical_distance(big.point, small.point)
    # small circle fits in big circle
    return dist + small.radius < big.radius 
end
function _contains(cap::SphericalCap, point::UnitSphericalPoint)
    spherical_distance(cap.point, point) <= cap.radius
end

#Comment by asinghvi: this could be transformed to GO.union
function _merge(x::SphericalCap, y::SphericalCap)

    d = spherical_distance(x.point, y.point)
    newradius = (x.radius + y.radius + d) / 2
    if newradius < x.radius
        #x contains y
        x
    elseif newradius < y.radius
        #y contains x
        y
    else
        excenter = 0.5 * (1 - (x.radius - y.radius) / d)
        newcenter = slerp(x.point, y.point, excenter)
        SphericalCap(newcenter, newradius)
    end
end

function circumcenter_on_unit_sphere(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    LinearAlgebra.normalize(a × b + b × c + c × a)
end

"Get the circumcenter of the triangle (a, b, c) on the unit sphere.  Returns a normalized 3-vector."
function SphericalCap(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    circumcenter = circumcenter_on_unit_sphere(a, b, c)
    circumradius = spherical_distance(a, circumcenter)
    return SphericalCap(circumcenter, circumradius)
end

function _is_ccw_unit_sphere(v_0::S, v_c::S, v_i::S) where S <: UnitSphericalPoint
    # checks if the smaller interior angle for the great circles connecting u-v and v-w is CCW
    return(LinearAlgebra.dot(LinearAlgebra.cross(v_c - v_0,v_i - v_c), v_i) < 0)
end

function angle_between(a::S, b::S, c::S) where S <: UnitSphericalPoint
    ab = b - a
    bc = c - b
    norm_dot = (ab ⋅ bc) / (LinearAlgebra.norm(ab) * LinearAlgebra.norm(bc))
    angle =  acos(clamp(norm_dot, -1.0, 1.0))
    if _is_ccw_unit_sphere(a, b, c)
        return angle
    else
        return 2π - angle
    end
end
