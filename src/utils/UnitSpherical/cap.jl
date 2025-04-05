#=
# Spherical caps
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
# TODO: add implementations to merge two spherical caps
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
