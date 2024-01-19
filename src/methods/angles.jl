# # Angles
export angles

#=
## What is angles?

Angles are the angles formed by a given geometries line segments, if it has line segments.

To provide an example, consider this rectangle:
```@example angles
using GeometryOps
using GeometryOps.GeometryBasics
using Makie

rect = Polygon([Point(0,0), Point(0,1), Point(1,1), Point(1,0), Point(0, 0)])
f, a, p = poly(rect; axis = (; aspect = DataAspect()))
```
This is clearly a rectangle, with angles of 90 degrees.
```@example angles
lines!(a, rect; color = 1:length(coordinates(rect))+1)
f
angles(area)
```

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a
wrapper method that dispatches to the correct implementation based on the
geometry trait. This is also used in the implementation, since it's a lot less
work!
=#

"""
    angles(geom, ::Type{T} = Float64)

Returns the angles of a geometry or collection of geometries. 
This is computed differently for different geometries:

    - The angles of a point is an empty vector.
    - The angles of a single line segment is an empty vector.
    - The angles of a linestring or linearring is a vector of angles formed by the curve.
    - The angles of a polygin is a vector of vectors of angles formed by each ring.
    - The angles of a multi-geometry collection is a vector of the angles of each of the
        sub-geometries as defined above.

Result will be a Vector, or nested set of vectors, of type T where an optional argument with
a default value of Float64.
"""
angles(geom, ::Type{T} = Float64) where T <: AbstractFloat =
    _angles(T, GI.trait(geom), geom)

# Points and single line segments have no angles
_angles(::Type{T}, ::Union{GI.PointTrait, GI.LineTrait}, geom) where T = T[]

#= The angles of a linestring are the angles formed by the line. If the first and last point
are not explicitly repeated, the geom is not considered closed. The angles should all be on
one side of the line, but a particular side is not guaranteed by this function. =#
function _angles(::Type{T}, ::Union{GI.LineStringTrait}, geom) where T
    npoints = GI.npoint(geom)
    first_last_equal = equals(GI.getpoint(geom, 1), GI.getpoint(geom, npoints))
    angle_list = Vector{T}(undef, npoints - (first_last_equal ? 1 : 2))
    _find_angles!(
        T, angle_list, geom;
        offset = first_last_equal, close_geom = false,
    )
    return angle_list
end

#= The angles of a linearring are the angles within the closed line and include the angles
formed by connecting the first and last points of the curve. =#
function _angles(::Type{T}, ::GI.LinearRingTrait, geom) where T
    npoints = GI.npoint(geom)
    first_last_equal = equals(GI.getpoint(geom, 1), GI.getpoint(geom, npoints))
    angle_list = Vector{T}(undef, npoints - (first_last_equal ? 1 : 0))
    _find_angles!(
        T, angle_list, geom;
        offset = true, close_geom = !first_last_equal,
    )
    return angle_list
end

#= The angles of a polygon is a vector of lists angles of its rings. Note that 
this means that the angles of any holes will be the exterior angles of the holes
outside of the geometry, rather than the interior angles=#
_angles(::Type{T}, ::GI.PolygonTrait, geom) where T =
    [_angles(T, GI.LinearRingTrait(), g) for g in GI.getring(geom)]

# Angles of a multi-geometry is simply a list of the angles of its sub-geometries.
_angles(::Type{T}, ::MultiGeomTrait, geom) where T = [angles(g, T) for g in GI.getgeom(geom)]

#=
Find angles of a curve and insert the values into the angle_list. If offset is true, then
save space for the angle at the first vertex, as the curve is closed, at the front of
angle_list. If close_geom is true, then despite the first and last point not being
explicitly repeated, the curve is closed and the angle of the last point should be added to
angle_list.
=#
function _find_angles!(::Type{T}, angle_list, geom; offset, close_geom) where T
    local p1, prev_p1_diff, p2_p1_diff
    local start_point, start_diff
    local extreem_idx, extreem_x, extreem_y
    i_offset = offset ? 1 : 0
    # Loop through the curve and find each of the angels
    for (i, p2) in enumerate(GI.getpoint(geom))
        xp2, yp2 = GI.x(p2), GI.y(p2)
        #= Find point with smallest x values (and smallest y in case of a tie) as this point
        is know to be convex. =#
        if i == 1 || (xp2 < extreem_x || (xp2 == extreem_x && yp2 < extreem_y))
            extreem_idx = i
            extreem_x, extreem_y = xp2, yp2
        end
        if i > 1
            p2_p1_diff = (xp2 - GI.x(p1), yp2 - GI.y(p1))
            if i == 2
                start_point = p1
                start_diff = p2_p1_diff
            else
                angle_list[i - 2 + i_offset] = _diffs_calc_angle(T, prev_p1_diff, p2_p1_diff)
            end
            prev_p1_diff = -1 .* p2_p1_diff
        end
        p1 = p2
    end
    # If the last point of geometry should be the same as the first, calculate closing angle
    if close_geom
        p2_p1_diff = (GI.x(start_point) - GI.x(p1), GI.y(start_point) - GI.y(p1))
        angle_list[end] = _diffs_calc_angle(T, prev_p1_diff, p2_p1_diff)
        prev_p1_diff = -1 .* p2_p1_diff
    end
    # If needed, calculate first angle corresponding to the first point 
    if offset
        angle_list[1] = _diffs_calc_angle(T, prev_p1_diff, start_diff)
    end
    #= Make sure that all of the angles are on the same side of the line and inside of the
    closed ring if the input geometry is closed. =#
    convex_sgn = sign(angle_list[extreem_idx])
    for i in eachindex(angle_list)
        idx_sgn = sign(angle_list[i])
        if idx_sgn == -1
            angle_list[i] = abs(angle_list[i])
        end
        if idx_sgn != convex_sgn
            angle_list[i] = 360 - angle_list[i]
        end
    end
    return
end

#=
Calculate the angle between two vectors defined by the previous and current Δx and Δys.
Angle will have a sign corresponding to the sign of the cross product between the two
vectors. All angles of one sign in a given geometry are convex, while those of the other
sign are concave. However, the sign corresponding to each of these can vary based on
geometry and thus you must compare to an angle that is know to be convex or concave.
=#
function _diffs_calc_angle(::Type{T}, (Δx_prev, Δy_prev), (Δx_curr, Δy_curr)) where T
    cross_prod = Δx_prev * Δy_curr - Δy_prev * Δx_curr
    dot_prod = Δx_prev * Δx_curr + Δy_prev * Δy_curr
    prev_mag = max(sqrt(Δx_prev^2 + Δy_prev^2), eps(T))
    curr_mag = max(sqrt(Δx_curr^2 + Δy_curr^2), eps(T))
    val = clamp(dot_prod / (prev_mag * curr_mag), -one(T), one(T))
    angle = real(acos(val) * 180 / π)
    return angle * (cross_prod < 0 ? -1 : 1)
end