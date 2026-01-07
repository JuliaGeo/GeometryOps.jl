using Makie
import GeometryOps as GO
import GeometryOps.UnitSpherical
import LinearAlgebra
import Makie.GeometryBasics
import StaticArrays


# Taken from Rotations.jl
function _angle_axis_rotation(theta, axis, v)
    # Using Rodrigues formula on an AngleAxis parametrization (assume unit axis length) to do the rotation
    # (implementation from: https://ceres-solver.googlesource.com/ceres-solver/+/1.10.0/include/ceres/rotation.h)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    w = axis ./ LinearAlgebra.norm(axis)
    st, ct = sincos(theta)
    w_cross_pt = LinearAlgebra.cross(w, v)
    m = LinearAlgebra.dot(v, w) * (one(w_cross_pt[1]) - ct)
    T = Base.promote_op(*, Base.promote_type(typeof(theta), eltype(axis)), eltype(v))
    return StaticArrays.similar_type(v,T)(v[1] * ct + w_cross_pt[1] * st + w[1] * m,
                                v[2] * ct + w_cross_pt[2] * st + w[2] * m,
                                v[3] * ct + w_cross_pt[3] * st + w[3] * m)
end

# Find a unit vector perpendicular to the given unit vector
function _perpendicular_vector(v)
    # Choose a reference vector that's not parallel to v
    if abs(v[3]) < 0.9
        ref = Point3d(0, 0, 1)
    else
        ref = Point3d(1, 0, 0)
    end
    return LinearAlgebra.normalize(LinearAlgebra.cross(v, ref))
end

# Generate a point on the boundary of a spherical cap
# theta: angle around the cap boundary (0 to 2π)
# phi: angular distance from cap center (0 to cap.radius)
# cap: the SphericalCap
function _spherical_cap_point(theta, phi, cap)
    if phi == 0
        return cap.point
    end

    # Get a perpendicular direction to the cap center
    perp = _perpendicular_vector(cap.point)

    # Rotate the perpendicular vector around cap.point by angle theta
    direction = _angle_axis_rotation(theta, cap.point, perp)
    direction = LinearAlgebra.normalize(direction)

    # Move from cap.point along this direction by angular distance phi
    # Using spherical geometry: point = cos(phi)*center + sin(phi)*direction
    result = cos(phi) * cap.point + sin(phi) * direction
    return LinearAlgebra.normalize(result)
end

function Makie.convert_arguments(::Type{Makie.Mesh}, cap::UnitSpherical.SphericalCap)
    N = 40
    rmin = 0.1
    rect = GeometryBasics.Tesselation(Rect2d(0, 0, 2pi, cap.radius), (N, max(2, ceil(Int, cap.radius / rmin))))
    faces = GeometryBasics.decompose(Makie.GLTriangleFace, rect)
    points = GeometryBasics.coordinates(rect)
    # Remove the first N points and set the first point to (0,0)
    # such that the first point is the center of the cap.
    points = points[N:end]
    points[1] = Makie.Point2d(0)
    faces = map(faces) do f
        Makie.GLTriangleFace(
            f[1] <= N ? Makie.GLIndex(1) : f[1] - N + 1,
            f[2] <= N ? Makie.GLIndex(1) : f[2] - N + 1,
            f[3] <= N ? Makie.GLIndex(1) : f[3] - N + 1,
        )
    end
    # Convert back to unit-spherical space.
    three_d_points = map(points) do (theta, phi)
        _spherical_cap_point(theta, phi, cap)
    end
    # Convert back to geographic space.
    # This is so that `zlevel` works right.
    return (GeometryBasics.normal_mesh(UnitSpherical.GeographicFromUnitSphere().(three_d_points) .|> Point2d, faces),)
    # return (GeometryBasics.normal_mesh(three_d_points, faces),)
end
Makie.convert_arguments(::Type{Makie.Poly}, cap::UnitSpherical.SphericalCap) = Makie.convert_arguments(Makie.Mesh, cap)
Makie.plottype(::Type{<:UnitSpherical.SphericalCap}) = Makie.Poly

function Makie.convert_arguments(::Makie.PointBased, cap::UnitSpherical.SphericalCap)
    N = 40
    points = GeometryBasics.Point2d.(LinRange(0, 2pi, N), cap.radius)

    three_d_points = map(points) do (theta, phi)
        _spherical_cap_point(theta, phi, cap)
    end
    # Convert back to geographic space.
    # This is so that `zlevel` works right.
    return (UnitSpherical.GeographicFromUnitSphere().(three_d_points) .|> Point2d,)
end