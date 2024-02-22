# # Area and signed area

export area, signed_area, coverage

#=
## What is area? What is signed area? What is coverage?

Area is the amount of space occupied by a two-dimensional figure. It is always
a positive value. Signed area is simply the integral over the exterior path of
a polygon, minus the sum of integrals over its interior holes. It is signed such
that a clockwise path has a positive area, and a counterclockwise path has a
negative area. The area is the absolute value of the signed area.

To provide an example, consider this rectangle:
```@example rect
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

rect = GI.Polygon([[(0,0), (0,1), (1,1), (1,0), (0, 0)]])
f, a, p = poly(collect(GI.getpoint(rect)); axis = (; aspect = DataAspect()))
```
This is clearly a rectangle, etc.  But now let's look at how the points look:
```@example rect
lines!(
    collect(GI.getpoint(rect));
    color = 1:GI.npoint(rect), linewidth = 10.0)
f
```
The points are ordered in a counterclockwise fashion, which means that the signed area
is negative.  If we reverse the order of the points, we get a postive area.
```@example rect
GO.signed_area(rect)  # -1.0
```

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a
wrapper method that dispatches to the correct implementation based on the
geometry trait. This is also used in the implementation, since it's a lot less
work!

Note that area (and signed area) are zero for all points and curves, even
if the curves are closed like with a linear ring. Also note that signed area
really only makes sense for polygons, given with a multipolygon can have several
polygons each with a different orientation and thus the absolute value of the
signed area might not be the area. This is why signed area is only implemented
for polygons.
=#

const _AREA_TARGETS = Union{GI.PolygonTrait,GI.AbstractCurveTrait,GI.MultiPointTrait,GI.PointTrait}

const UNKNOWN, NORTH, EAST, SOUTH, WEST = 0:4
"""
    area(geom, ::Type{T} = Float64)::T

Returns the area of a geometry or collection of geometries. 
This is computed slightly differently for different geometries:

    - The area of a point/multipoint is always zero.
    - The area of a curve/multicurve is always zero.
    - The area of a polygon is the absolute value of the signed area.
    - The area multi-polygon is the sum of the areas of all of the sub-polygons.
    - The area of a geometry collection, feature collection of array/iterable 
        is the sum of the areas of all of the sub-geometries. 

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function area(geom, ::Type{T} = Float64; threaded=false) where T <: AbstractFloat
    applyreduce(+, _AREA_TARGETS, geom; threaded, init=zero(T)) do g
        _area(T, GI.trait(g), g)
    end
end


"""
    signed_area(geom, ::Type{T} = Float64)::T

Returns the signed area of a single geometry, based on winding order. 
This is computed slighly differently for different geometries:

    - The signed area of a point is always zero.
    - The signed area of a curve is always zero.
    - The signed area of a polygon is computed with the shoelace formula and is
    positive if the polygon coordinates wind clockwise and negative if
    counterclockwise.
    - You cannot compute the signed area of a multipolygon as it doesn't have a
    meaning as each sub-polygon could have a different winding order.

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
signed_area(geom, ::Type{T} = Float64) where T <: AbstractFloat =
    _signed_area(T, GI.trait(geom), geom)

# Points, MultiPoints, Curves, MultiCurves
_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

_signed_area(::Type{T}, ::GI.AbstractGeometryTrait, geom) where T = zero(T)

# Polygons
_area(::Type{T}, trait::GI.PolygonTrait, poly) where T =
    abs(_signed_area(T, trait, poly))

function _signed_area(::Type{T}, ::GI.PolygonTrait, poly) where T
    GI.isempty(poly) && return zero(T)
    s_area = _signed_area(T, GI.getexterior(poly))
    area = abs(s_area)
    area == 0 && return area
    # Remove hole areas from total
    for hole in GI.gethole(poly)
        area -= abs(_signed_area(T, hole))
    end
    # Winding of exterior ring determines sign
    return area * sign(s_area)
end

#=
Helper function:

Calculates the signed area of a given curve. This is equivalent to integrating
to find the area under the curve. Even if curve isn't explicitly closed by
repeating the first point at the end of the coordinates, curve is still assumed
to be closed.
=#
function _signed_area(::Type{T}, geom) where T
    area = zero(T)
    np = GI.npoint(geom)
    np == 0 && return area

    first = true
    local pfirst, p1
    # Integrate the area under the curve
    for p2 in GI.getpoint(geom)
        # Skip the first and do it later 
        # This lets us work within one iteration over geom, 
        # which means on C call when using points from external libraries.
        if first
            p1 = pfirst = p2
            first = false
            continue
        end
        # Accumulate the area into `area`
        area += _area_component(p1, p2)
        p1 = p2
    end
    # Complete the last edge.
    # If the first and last where the same this will be zero
    p2 = pfirst
    area += _area_component(p1, p2)
    return T(area / 2)
end

# One term of the shoelace area formula
_area_component(p1, p2) = GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)

"""
    calc_floe_area_in_cell(coords, xmin, xmax, ymin, ymax)

Calculate area of intersection between given polygon and grid cell defined by
its minimum and maximum x and y-coordinates. 
Inputs:
    coords  <PolyVec> polygon coordinates
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
Output:
    Area of intersection between polygon and grid cell
Warning:
    Assumes polygon has clockwise winding order
"""
coverage(geom, xmin, xmax, ymin, ymax, ::Type{T} = Float64) where T =
    _coverage(T, GI.trait(geom), geom, xmin, xmax, ymin, ymax)

# Points, MultiPoints, Curves, MultiCurves
_coverage(::Type{T}, ::GI.AbstractGeometryTrait, geom, cell_extremes...) where T = zero(T)

# Polygons
function _coverage(::Type{T}, ::GI.PolygonTrait, poly, cell_extremes...) where T
    GI.isempty(poly) && return zero(T)
    cov_area = _coverage(T, GI.getexterior(poly), cell_extremes...)
    cov_area == 0 && return cov_area
    # Remove hole areas from total
    for hole in GI.gethole(poly)
        cov_area -= _coverage(T, hole, cell_extremes...)
    end
    # Winding of exterior ring determines sign
    return cov_area
end

#=
Helper functions:

...
=#
function _coverage(::Type{T}, ring, xmin, xmax, ymin, ymax) where T
    cov_area = zero(T)
    unmatched_out_wall, unmatched_out_point = UNKNOWN, (zero(T), zero(T))
    unmatched_in_wall, unmatched_in_point = unmatched_out_wall, unmatched_out_point
    # Loop over edges of polygon
    ring_cw = isclockwise(ring)
    p1 = _tuple_point(GI.getpoint(ring, ring_cw ? GI.npoint(ring) : 1))
    for p in (ring_cw ? GI.getpoint(ring) : reverse(GI.getpoint(ring)))
        p2 = _tuple_point(p)
        # Determine if edge points are within the cell
        p1_in_cell = _point_in_cell(p1, xmin, xmax, ymin, ymax)
        p2_in_cell = _point_in_cell(p2, xmin, xmax, ymin, ymax)
        # If entire line segment is inside cell
        if p1_in_cell && p2_in_cell
            cov_area += _area_component(p1, p2)
            p1 = p2
            continue
        end
        # If edge passes outside of rectangle, determine which edge segments are added
        inter1, inter2 = _line_intersect_cell(T, p1, p2, xmin, xmax, ymin, ymax)
        # Endpoints of segment within the cell and wall they are on if known
        (in_wall, in_point), (out_wall, out_point) =
            if p1_in_cell
                ((UNKNOWN, p1), inter1)
            elseif p2_in_cell
                (inter1, (UNKNOWN, p2))
            else
                distance(inter1[2], p1) < distance(inter2[2], p1) ?
                    (inter1, inter2) : (inter2, inter1)
            end
        # Add edge component
        cov_area += _area_component(in_point, out_point)
        # Connect intersection along cell walls
        if in_wall != UNKNOWN
            if unmatched_out_wall == UNKNOWN
                unmatched_in_point = in_point
                unmatched_in_wall = in_wall
            else
                cov_area += connect_edges(T, unmatched_out_point, in_point,
                    unmatched_out_wall, in_wall, xmin, xmax, ymin, ymax)
                unmatched_out_wall = out_wall
            end
        end
        if out_wall != UNKNOWN
            unmatched_out_wall, unmatched_out_point = out_wall, out_point
        end
        p1 = p2
    end
    # if unmatched in-point at beginning, close polygon with last out point
    if unmatched_in_wall != UNKNOWN
        cov_area += connect_edges(T, unmatched_out_point, unmatched_in_point,
            unmatched_out_wall, unmatched_in_wall, xmin, xmax, ymin, ymax)
    end
    cov_area = abs(cov_area) / 2
    #  if grid cell is within polygon then the area is grid cell area
    if cov_area == 0 && _point_filled_curve_orientation((xmin, ymin), ring;
        in = true, on = true, out = false)
        cov_area = abs((xmax - xmin) * (ymax - ymin))
    end
    return cov_area
end

_point_in_cell((x, y), xmin, xmax, ymin, ymax) = xmin <= x <= xmax && ymin <= y <= ymax

_between(b, c, a) = a ≤ b ≤ c || c ≤ b ≤ a 

function _line_intersect_cell(::Type{T}, (x1, y1), (x2, y2), xmin, xmax, ymin, ymax) where T
    Δx, Δy = x2 - x1, y2 - y1
    inter1 = (UNKNOWN, (zero(T), zero(T)))
    inter2 = inter1
    if Δx == 0
        if xmin ≤ x1 ≤ xmax
            inter1 = _between(ymax, y1, y2) ? (NORTH, (x1, ymax)) : inter1
            inter2 = _between(ymin, y1, y2) ? (SOUTH, (x1, ymin)) : inter2
        end
    elseif Δy == 0
        if ymin ≤ y1 ≤ ymax
            inter1 = _between(xmax, x1, x2) ? (EAST, (xmax, y1)) : inter1
            inter2 = _between(xmin, x1, x2) ? (WEST, (xmin, y1)) : inter2
        end
    else
        m = Δy / Δx
        b = y1 - m * x1
        # Calculate potential intersections
        xn = (ymax - b) / m
        if xmin ≤ xn ≤ xmax && _between(xn, x1, x2) && _between(ymax, y1, y2)
            inter1 = (NORTH, (xn, ymax))
        end
        xs = (ymin - b) / m
        if xmin ≤ xs ≤ xmax && _between(xs, x1, x2) && _between(ymin, y1, y2)
            new_intr = (SOUTH, (xs, ymin))
            (inter1[1] == UNKNOWN) ? (inter1 = new_intr) : (inter2 = new_intr)
        end
        ye =  m * xmax + b
        if ymin ≤ ye ≤ ymax && _between(ye, y1, y2) && _between(xmax, x1, x2)
            new_intr = (EAST, (xmax, ye))
            (inter1[1] == UNKNOWN) ? (inter1 = new_intr) : (inter2 = new_intr)
        end
        yw = m * xmin + b
        if ymin ≤ yw ≤ ymax && _between(yw, y1, y2) && _between(xmin, x1, x2)
            new_intr = (WEST, (xmin, yw))
            (inter1[1] == UNKNOWN) ? (inter1 = new_intr) : (inter2 = new_intr)
        end
    end
    if inter1[1] == UNKNOWN
        inter1, inter2 = inter2, inter1
    end
    return inter1, inter2
end

#=
    connect_edges(x1, y1, x2, y2, wall1, wall2, xmin, xmax, ymin, ymax)

    Area component of shoelace formula coming from the distance between point 1
    and point 2 along grid cell walls
=#
function connect_edges(::Type{T}, p1, p2, wall1, wall2, xmin, xmax, ymin, ymax) where {T}
    connect_area = zero(T)
    if wall1 == wall2 && _is_clockwise_from(p1, p2, wall1)
        connect_area += _area_component(p1, p2)
    else
        # From the point to the corner of wall 1
        connect_area += _partial_edge_out(p1, xmin, xmax, ymin, ymax, wall1)
        # Any intermediate walls (full length)
        next_wall, last_wall = wall1 + 1, wall2 - 1
        if wall2 > wall1
            for wall in next_wall:last_wall
                connect_area += _full_edge_area(xmin, xmax, ymin, ymax, wall)
            end
        else
            for wall in Iterators.flatten((next_wall:WEST, NORTH:last_wall))
                connect_area += _full_edge_area(xmin, xmax, ymin, ymax, wall)
            end
        end
        # From the corner of wall 2 to the point
        connect_area += _partial_edge_in(p2, xmin, xmax, ymin, ymax, wall2)
    end
    return connect_area
end

_is_clockwise_from((x1, y1), (x2, y2), wall) = (wall == NORTH && x2 > x1) ||
    (wall == EAST && y2 < y1) || (wall == SOUTH && x2 < x1) || (wall == WEST && y2 > y1)

_full_edge_area(xmin, xmax, ymin, ymax, wall) = if wall == NORTH
        ymax * (xmin - xmax)
    elseif wall == EAST
        xmax * (ymin - ymax)
    elseif wall == SOUTH
        ymin * (xmax - xmin)
    else
        xmin * (ymax - ymin)
    end

function _partial_edge_in((x2, y2), xmin, xmax, ymin, ymax, wall)
    x_wall = (wall == NORTH || wall == WEST) ? xmin : xmax
    y_wall = (wall == NORTH || wall == EAST) ? ymax : ymin
    return x_wall * y2 - x2 * y_wall
end

function _partial_edge_out((x1, y1), xmin, xmax, ymin, ymax, wall)
    x_wall = (wall == NORTH || wall == EAST) ? xmax : xmin
    y_wall = (wall == NORTH || wall == WEST) ? ymax : ymin
    return x1 * y_wall - x_wall * y1
end
