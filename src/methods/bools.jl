# # Boolean conditions

export isclockwise, isconcave
export point_on_line, point_in_polygon, point_in_ring
export line_on_line, line_in_polygon, polygon_in_polygon

# These are all adapted from Turf.jl.

# The may not necessarily be what want in the end but work for now!

"""
    isclockwise(line::Union{LineString, Vector{Position}})::Bool

Take a ring and return true or false whether or not the ring is clockwise or counter-clockwise.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line = GI.LineString([(0, 0), (1, 1), (1, 0), (0, 0)])
GO.isclockwise(line)

# output
true
```
"""
isclockwise(geom)::Bool = isclockwise(GI.trait(geom), geom)
function isclockwise(::AbstractCurveTrait, line)::Bool
    sum = 0.0
    prev = GI.getpoint(line, 1)
    for p in GI.getpoint(line)
        # sum will be zero for the first point as x is subtracted from itself
        sum += (GI.x(p) - GI.x(prev)) * (GI.y(p) + GI.y(prev))
        prev = p
    end

    return sum > 0.0
end

"""
    isconcave(poly::Polygon)::Bool

Take a polygon and return true or false as to whether it is concave or not.

## Examples
```jldoctest
import GeoInterface as GI, GeometryOps as GO

poly = GI.Polygon([[(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]])
GO.isconcave(poly)

# output
false
```
"""
function isconcave(poly)::Bool
    sign = false

    exterior = GI.getexterior(poly)

    # FIXME handle not closed polygons
    GI.npoint(exterior) <= 4 && return false
    n = GI.npoint(exterior) - 1

    for i in 1:n
        j = ((i + 1) % n) === 0 ? 1 : (i + 1) % n
        m = ((i + 2) % n) === 0 ? 1 : (i + 2) % n

        pti = GI.getpoint(exterior, i)
        ptj = GI.getpoint(exterior, j)
        ptm = GI.getpoint(exterior, m)

        dx1 = GI.x(ptm) - GI.x(ptj)
        dy1 = GI.y(ptm) - GI.y(ptj)
        dx2 = GI.x(pti) - GI.x(ptj)
        dy2 = GI.y(pti) - GI.y(ptj)

        cross = (dx1 * dy2) - (dy1 * dx2)

        if i === 0
            sign = cross > 0
        elseif sign !== (cross > 0)
            return true
        end
    end

    return false
end

equals(geo1, geo2) = _equals(trait(geo1), geo1, trait(geo2), geo2)

_equals(::T, geo1, ::T, geo2) where T = error("Cant compare $T yet")
function _equals(::T, p1, ::T, p2) where {T<:PointTrait}
    GI.ncoord(p1) == GI.ncoord(p2) || return false
    GI.x(p1) == GI.x(p2) || return false
    GI.y(p1) == GI.y(p2) || return false
    if GI.is3d(p1)
        GI.z(p1) == GI.z(p2) || return false 
    end
    return true
end
function _equals(::T, l1, ::T, l2) where {T<:AbstractCurveTrait}
    # Check line lengths match
    GI.npoint(l1) == GI.npoint(l2) || return false

    # Then check all points are the same
    for (p1, p2) in zip(GI.getpoint(l1), GI.getpoint(l2))
        equals(p1, p2) || return false
    end
    return true
end
_equals(t1, geo1, t2, geo2) = false

# """
#     isparallel(line1::LineString, line2::LineString)::Bool

# Return `true` if each segment of `line1` is parallel to the correspondent segment of `line2`

# ## Examples
# ```julia
# import GeoInterface as GI, GeometryOps as GO
# julia> line1 = GI.LineString([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)])
# GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.170356, 45.477985), (9.164434, 45.482551), (9.166644, 45.484003)], nothing, nothing)

# julia> line2 = GI.LineString([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)])
# GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(9.169356, 45.477985), (9.163434, 45.482551), (9.165644, 45.484003)], nothing, nothing)

# julia> 
# GO.isparallel(line1, line2)
# true
# ```
# """
# function isparallel(line1, line2)::Bool
#     seg1 = linesegment(line1)
#     seg2 = linesegment(line2)

#     for i in eachindex(seg1)
#         coors2 = nothing
#         coors1 = seg1[i]
#         coors2 = seg2[i]
#         _isparallel(coors1, coors2) == false && return false
#     end
#     return true
# end

# @inline function _isparallel(p1, p2)
#     slope1 = bearing_to_azimuth(rhumb_bearing(GI.x(p1), GI.x(p2)))
#     slope2 = bearing_to_azimuth(rhumb_bearing(GI.y(p1), GI.y(p2)))

#     return slope1 === slope2
# end


"""
    point_on_line(point::Point, line::LineString; ignore_end_vertices::Bool=false)::Bool

Return true if a point is on a line. Accept a optional parameter to ignore the
start and end vertices of the linestring.

## Examples

```jldoctest
import GeoInterface as GI, GeometryOps as GO

point = (1, 1)
line = GI.LineString([(0, 0), (3, 3), (4, 4)])
GO.point_on_line(point, line)

# output
true
```
"""
function point_on_line(point, line; ignore_end_vertices::Bool=false)::Bool
    line_points = tuple_points(line)
    n = length(line_points)

    exclude_boundary = :none
    for i in 1:n - 1
        if ignore_end_vertices
            if i === 1 
                exclude_boundary = :start
            elseif i === n - 2 
                exclude_boundary = :end
            elseif (i === 1 && i + 1 === n - 1) 
                exclude_boundary = :both
            end
        end
        if point_on_segment(point, (line_points[i], line_points[i + 1]); exclude_boundary) 
            return true
        end
    end
    return false
end

function point_on_segment(point, (start, stop); exclude_boundary::Symbol=:none)::Bool
    x, y = GI.x(point), GI.y(point)
    x1, y1 = GI.x(start), GI.y(start)
    x2, y2 = GI.x(stop), GI.y(stop)

    dxc = x - x1
    dyc = y - y1
    dx1 = x2 - x1
    dy1 = y2 - y1

    # TODO use better predicate for crossing here
    cross = dxc * dy1 - dyc * dx1
    cross != 0 && return false

    # Will constprop optimise these away?
    if exclude_boundary === :none
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x <= x2 : x2 <= x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y <= y2 : y2 <= y && y <= y1
    elseif exclude_boundary === :start
        if abs(dx1) >= abs(dy1)
             return dx1 > 0 ? x1 < x && x <= x2 : x2 <= x && x < x1
        end
        return dy1 > 0 ? y1 < y && y <= y2 : y2 <= y && y < y1
    elseif exclude_boundary === :end
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x < x2 : x2 < x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y < y2 : y2 < y && y <= y1
    elseif exclude_boundary === :both
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 < x && x < x2 : x2 < x && x < x1
        end
        return dy1 > 0 ? y1 < y && y < y2 : y2 < y && y < y1
    end
    return false
end

"""
    point_in_polygon(point::Point, polygon::Union{Polygon, MultiPolygon}, ignore_boundary::Bool=false)::Bool

Take a Point and a Polygon and determine if the point
resides inside the polygon. The polygon can be convex or concave. The function accounts for holes.

## Examples

```jldoctest
import GeoInterface as GI, GeometryOps as GO

point = (-77.0, 44.0)
poly = GI.Polygon([[(-81, 41), (-81, 47), (-72, 47), (-72, 41), (-81, 41)]])
GO.point_in_polygon(point, poly)

# output
true
```
"""
point_in_polygon(point, polygon; kw...)::Bool =
    point_in_polygon(GI.trait(point), point, GI.trait(polygon), polygon; kw...)
function point_in_polygon(
    ::PointTrait, point, 
    ::PolygonTrait, poly; 
    ignore_boundary::Bool=false,
    check_extent::Bool=false,
)::Bool
    # Cheaply check that the point is inside the polygon extent
    if check_extent
        point_in_extent(point, GI.extent(poly)) || return false
    end

    # Then check the point is inside the exterior ring
    point_in_polygon(point, GI.getexterior(poly); ignore_boundary, check_extent=false) || return false

    # Finally make sure the point is not in any of the holes,
    # flipping the boundary condition
    for ring in GI.gethole(poly)
        point_in_polygon(point, ring; ignore_boundary=!ignore_boundary) && return false
    end
    return true
end
function point_in_polygon(
    ::PointTrait, pt, 
    ::Union{LineStringTrait,LinearRingTrait}, ring; 
    ignore_boundary::Bool=false,
    check_extent::Bool=false,
)::Bool
    # Cheaply check that the point is inside the ring extent
    if check_extent
        point_in_extent(point, GI.extent(ring)) || return false
    end

    # Then check the point is inside the ring
    inside = false
    n = GI.npoint(ring)
    p_start = GI.getpoint(ring, 1)
    p_end = GI.getpoint(ring, n)

    # Handle closed on non-closed rings
    l = if GI.x(p_start) == GI.x(p_end) && GI.y(p_start) == GI.y(p_end) 
        l = n - 1
    else
        n
    end

    # Loop over all points in the ring
    for i in 1:l - 1
        j = i + 1

        p_i = GI.getpoint(ring, i)
        p_j = GI.getpoint(ring, j)
        xi = GI.x(p_i)
        yi = GI.y(p_i)
        xj = GI.x(p_j)
        yj = GI.y(p_j)

        on_boundary = (GI.y(pt) * (xi - xj) + yi * (xj - GI.x(pt)) + yj * (GI.x(pt) - xi) == 0) &&
            ((xi - GI.x(pt)) * (xj - GI.x(pt)) <= 0) && ((yi - GI.y(pt)) * (yj - GI.y(pt)) <= 0)

        on_boundary && return !ignore_boundary

        intersects = ((yi > GI.y(pt)) !== (yj > GI.y(pt))) && 
            (GI.x(pt) < (xj - xi) * (GI.y(pt) - yi) / (yj - yi) + xi)

        if intersects 
            inside = !inside
        end
    end

    return inside
end

function point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y1) = extent.X, extent.Y
    return x1 <= GI.x(p) && y1 <= GI.y(p) && x2 >= GI.x(p) && y2 >= GI.y(p)
end

function line_in_polygon(line, polygon)
    out = false

    Extents.intersects(GI.extent(polygon), GI.extent(line)) || return false

    p1 = GI.getpoint(line, 1)

    for i in 1:GI.npoint(line)
        p2 = GI.getpoint(line, i)
        mid = (GI.x(p1) + GI.x(p2)) / 2, (GI.y(p1) + GI.y(p2)) / 2

        # FIXME mid point in the polygon? what is that testing?
        if point_in_polygon(mid, poly; ignore_boundary=true)
            out = true
            break
        end
        p1 = p2
    end
    return out
end

line_on_line(line1, line2) = line_on_line(trait(line1), line1, trait(line2), line2)
function line_on_line(t1::GI.AbstractCurveTrait, line1, t2::AbstractCurveTrait, line2)
    for p in GI.getpoint(line1)
        # FIXME: all points being on the line doesn't
        # actually mean the whole line is on the line...
        point_on_line(p, line2) || return false
    end
    return true
end

line_in_polygon(line, poly) = line_in_polygon(trait(line), line, trait(poly), poly)
function line_in_polygon(
    ::LineStringTrait, line, 
    ::Union{AbstractPolygonTrait,LinearRingTrait}, poly
)
    Extents.intersects(GI.extent(poly), GI.extent(line)) && return false

    inside = false
    for i in 1:GI.npoint(line) - 1
        p = GI.getpont(line, i)
        p2 = GI.getpont(line, i + 1)
        point_in_polygon(p, poly) || return false
        if !inside 
            inside = point_in_polygon(p, poly; ignore_boundary=true)
        end
        # FIXME This seems like a hack, we should check for intersections rather than midpoint??
        if !inside
            mid = ((GI.x(p) + GI.x(p2)) / 2, (GI.y(p) + GI.y(p2)) / 2)
            inside = point_in_polygon(Point(mid), poly; ignore_boundary=true)
        end
    end
    return inside
end

function polygon_in_polygon(poly1, poly2)
     # Check the extents intersect
     Extents.intersects(GI.extent(poly1), GI.extent(poly2)) || return false
     # Check all points in poly1 are in poly2
     for point in GI.getpoint(poly1)
         point_in_polygon(point, poly2) || return false
     end
     # Check no points in poly2 are in poly1
     for point in GI.getpoint(poly2)
         point_in_polygon(point, poly1; ignore_boundary=true) && return false
     end
     # Check poly1 does not intersect poly2
     intersects(poly1, poly2) && return false

     # poly1 must be in poly2
     return true
 end
