# # Area and signed area

export area, signed_area

#=
## What is area? What is signed area?

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
        area += GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)
        p1 = p2
    end
    # Complete the last edge.
    # If the first and last where the same this will be zero
    p2 = pfirst
    area += GI.x(p1) * GI.y(p2) - GI.y(p1) * GI.x(p2)
    return T(area / 2)
end

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
function calc_floe_area_in_cell(
    coords::PolyVec{FT},
    xmin, xmax, ymin, ymax,
) where {FT <: AbstractFloat}
    npoints = length(coords[1]) - 1
    sum = FT(0)
    # keep track of when polygon edges leave grid cell
    out_point::Union{Nothing, Tuple{FT, FT}} = nothing
    out_wall = 0
    #=
    keep track of when polygon edges enter grid cell if it doesn't have
    matching outpoint from when polygon left grid cell
    =#
    unmatched_in_point = out_point
    unmatched_in_wall = 0
    # keep track of wall intersections
    iwalls = [FT(Inf) for _ in 1:4]
    @views for i in 1:npoints
        edge_area = FT(0)
        x1, y1 = coords[1][i]
        x2, y2 = coords[1][i + 1]
        # check if points are withing grid cell
        p1_in_cell = point_in_cell(x1, y1, xmin, xmax, ymin, ymax)
        p2_in_cell = point_in_cell(x2, y2,xmin, xmax, ymin, ymax)
        if p1_in_cell && p2_in_cell  # line segment inside cell
            edge_area += x1 * y2 - x2 * y1
        else  # at least one point outside of cell (intersection!)
            line_intersect_cell!(x1, y1, x2, y2, xmin, xmax, ymin, ymax, iwalls)
            i1_idx = findfirst(!isinf, iwalls)
            if p1_in_cell  # line segment exits cell
                out_point = get_inter_point(iwalls, i1_idx, xmin, xmax, ymin, ymax)
                out_wall = i1_idx
                edge_area += x1 * out_point[2] - out_point[1] * y1
            elseif p2_in_cell  # line segment enters cell
                in_point = get_inter_point(iwalls, i1_idx, xmin, xmax, ymin, ymax)
                edge_area += in_point[1] * y2 - x2 * in_point[2]
                if isnothing(out_point)  # need to match with outpoint at end
                    unmatched_in_point = in_point
                    unmatched_in_wall = i1_idx
                else # add edge areas connecting last out_point to this in_point
                    edge_area += connect_edges(
                        out_point[1], out_point[2],
                        in_point[1], in_point[2],
                        out_wall, i1_idx,
                        xmin, xmax, ymin, ymax,
                    )
                    out_point = nothing
                end
            elseif !isnothing(i1_idx)  # line passes through cell (in and out)
                i1 = get_inter_point(iwalls, i1_idx, xmin, xmax, ymin, ymax)
                i2_idx = findnext(!isinf, iwalls, i1_idx + 1)
                i2 = get_inter_point(iwalls, i2_idx, xmin, xmax, ymin, ymax)
                # Determine which intersection goes with which point
                x1_to_i1_dist = sqrt((i1[1] - x1)^2 + (i1[2] - y1)^2)
                x1_to_i2_dist = sqrt((i2[1] - x1)^2 + (i2[2] - y1)^2)
                if x1_to_i1_dist > x1_to_i2_dist
                    i1, i2 = i2, i1
                    i1_idx, i2_idx = i2_idx, i1_idx
                end
                # Component from point i1 to i2
                edge_area += i1[1] * i2[2] - i2[1] * i1[2]
                if isnothing(out_point)  # need to match with outpoint at end 
                    unmatched_in_point = i1
                    unmatched_in_wall = i1_idx
                else  # connect last out to new in (i1)
                    edge_area += connect_edges(
                        out_point[1], out_point[2],
                        i1[1], i1[2],
                        out_wall, i1_idx,
                        xmin, xmax, ymin, ymax,
                    )
                end
                # keep track of polygon leaving cell
                out_point = i2
                out_wall = i2_idx
            end
        end
        sum += edge_area
    end
    # if unmatched in-point at beginning, close polygon with last out point
    if !isnothing(unmatched_in_point)
        sum += connect_edges(
            out_point[1], out_point[2],
            unmatched_in_point[1], unmatched_in_point[2],
            out_wall, unmatched_in_wall,
            xmin, xmax, ymin, ymax,
        )
    end
    sum = abs(sum) / 2
    #  if grid cell is within polygon then the area is grid cell area
    if sum == 0 && is_point_in_poly(xmin, ymin, coords)
        sum = abs((xmax - xmin) * (ymax - ymin))
    end
    return sum
end

"""
    point_in_cell(x, y, xmin, xmax, ymin, ymax)

Checks if a point (x, y) is within a grid cell defined by its maximum and
minimum corner values in the x and y-directions.
Inputs:
    x       <AbstractFloat> x-coordinate of point
    y       <AbstractFloat> y-coordinate of point
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
Outputs:
    True if (x, y) is within grid cell, false otherwise.
"""
point_in_cell(x, y, xmin, xmax, ymin, ymax) =
    xmin <= x <= xmax && ymin <= y <= ymax

"""
    which_cell_wall(x, y, xmin, xmax, ymin, ymax)

Determine which cell wall a point is on.
Inputs:
    x       <AbstractFloat> x-coordinate of point
    y       <AbstractFloat> y-coordinate of point
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
Outputs:
    '1' if point is on north wall, '2' if point is on east wall, '3' if point is
    on south wall, '4' if point is on west side, and '0' otherwise
"""
function which_cell_wall(x, y, xmin, xmax, ymin, ymax) 
    wall =
        if x == xmin
            4  # west
        elseif x == xmax
            2  # east
        elseif y == ymin
            3  # south
        elseif y == ymax
            1  # north
        else
            0  # not on a wall
        end
    return wall
end

"""
    line_intersect_cell!(x1, y1, x2, y2, xmin, xmax, ymin, ymax, iwalls)

Inputs:
    x1      <AbstractFloat> x-coordinate of first point
    y1      <AbstractFloat> y-coordinate of first point
    x2      <AbstractFloat> x-coordinate of second point
    y2      <AbstractFloat> y-coordinate of second point
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
    iwalls  <Vector{AbstractFloat}> vector for intersection points on north,
                east, south, and west walls - since at least one of coordinates
                is known by position on wall, only store unknown value
Outputs:
    Nothing. Updates iwalls with intersection values. Any indices with Inf
    values don't have an intersection point.
"""
function line_intersect_cell!(
    x1::FT, y1, x2, y2,
    xmin, xmax, ymin, ymax, iwalls,
) where {FT}
    Δx = x2 - x1
    Δy = y2 - y1
    if Δx == 0
        iwalls[2] = FT(Inf)
        iwalls[4] = FT(Inf)
        atol = min(xmax - xmin, Δy, ymax - ymin) * 10eps()
        if between(x1, xmin, xmax, atol) # line is vertical
            iwalls[1] = between(ymax, y1, y2, atol) ? x1 : FT(Inf)  # north
            iwalls[3] = between(ymin, y1, y2, atol) ? x1 : FT(Inf)  # south
        else
            iwalls[1] = FT(Inf)
            iwalls[3] = FT(Inf)
        end
    elseif Δy == 0
        iwalls[1] = FT(Inf)
        iwalls[3] = FT(Inf)
        atol = min(Δx, xmax - xmin, ymax - ymin) * 1e-8
        if between(y1, ymin, ymax, atol)  # line is horizontal
            iwalls[2] = between(xmax, x1, x2, atol) ? y1 : FT(Inf)  # east
            iwalls[4] = between(xmin, x1, x2, atol) ? y1 : FT(Inf)  # west
        else
            iwalls[2] = FT(Inf)
            iwalls[4] = FT(Inf)
        end
    else  # line is slanted
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
        # Calculate potential intersections
        x_north = (ymax - b) / m
        y_east = m * xmax + b
        x_south = (ymin - b) / m
        y_west = m * xmin + b
        # Determine if they really are intersections
        atol = min(Δx, xmax - xmin, Δy, ymax - ymin) * 1e-8
        iwalls[1] = between(x_north, x1, x2, atol) &&
            between(x_north, xmin, xmax, atol) &&
            between(ymax, y1, y2, atol) ?
            x_north :
            FT(Inf)
        iwalls[2] = between(xmax, x1, x2, atol) &&
            between(y_east, y1, y2, atol) &&
            between(y_east, ymin, ymax, atol) ?
            y_east :
            FT(Inf)
        iwalls[3] = between(x_south, x1, x2, atol) &&
            between(x_south, xmin, xmax, atol) &&
            between(ymin, y1, y2, atol) ?
            x_south :
            FT(Inf)
        iwalls[4] = between(xmin, x1, x2, atol) &&
            between(y_west, y1, y2, atol) &&
            between(y_west, ymin, ymax, atol) ?
            y_west :
            FT(Inf)
    end
    return
end

"""
    get_inter_point(iwalls, i, xmin, xmax, ymin, ymax)

Returns x and y of intersection point for wall i
Inputs:
    iwalls  <Vector{AbstractFloat}> vector for intersection points on north,
                east, south, and west walls - since at least one of coordinates
                is known by position on wall, only store unknown value
    i       <Int> wall index (1-4) of intersection point
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
Outputs:
    Intersection point with given wall    
"""
function get_inter_point(iwalls, i, xmin, xmax, ymin, ymax)
    if i == 1
        return iwalls[1], ymax
    elseif i == 2
        return xmax, iwalls[2]
    elseif i == 3
        return iwalls[3], ymin
    elseif i == 4
        return xmin, iwalls[4]
    end
end

"""
    is_clockwise_from(x1, y1, x2, y2, wall)

Checks if first point is clockwise from second point on the same wall
Inputs:
    x1      <Number> x-coordinate of first point
    y1      <Number> y-coordinate of first point
    x2      <Number> x-coordinate of second point
    y2      <Number> y-coordinate of second point
    wall    <Int> wall index (1-4)
"""
function is_clockwise_from(x1, y1, x2, y2, wall)
    if wall == 1  # north
        return x2 > x1
    elseif wall == 2  # east
        return y2 < y1
    elseif wall == 3  # south
        return x2 < x1
    elseif wall == 4  # west
        return y2 > y1
    end
end

"""
    full_edge_area(xmin, xmax, ymin, ymax, wall)

Calculate entire edge area component using shoelace formula
Inputs:
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
    wall    <Int> wall index (1-4)
Output:
    Entire edge area component using the shoelace formula
"""
function full_edge_area(xmin, xmax, ymin, ymax, wall)
    if wall == 1
        return ymax * (xmin - xmax)
    elseif wall == 2
        return xmax * (ymin - ymax)
    elseif wall == 3
        return ymin * (xmax - xmin)
    else
        return xmin * (ymax - ymin)
    end
end

"""
    partial_edge_in(x2, y2, xmin, xmax, ymin, ymax, wall)

Calculate area component from the corner of the wall to the given point using
the shoelace formula
Inputs:
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
    wall    <Int> wall index (1-4)
Output:
    Edge area component from the corner of the wall to the given point using
    the shoelace formula
"""
function partial_edge_in(x2, y2, xmin, xmax, ymin, ymax, wall)
    # from the corner to the point
    if wall == 1
        return xmin * y2 - x2 * ymax
    elseif wall == 2
        return xmax * y2 - x2 * ymax
    elseif wall == 3
        return xmax * y2 - x2 * ymin
    else
        return xmin * y2 - x2 * ymin
    end
end

"""
    partial_edge_out(x1, y1, xmin, xmax, ymin, ymax, wall)

Calculate area component from the point to the next corner of the wall using the
shoelace formula
Inputs:
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
    wall    <Int> wall index (1-4)
Output:
    Edge area component from the point to the next corner of the wall using the
    shoelace formula
"""
function partial_edge_out(x1, y1, xmin, xmax, ymin, ymax, wall)
    # from the corner to the point
    if wall == 1
        return x1 * ymax - xmax * y1
    elseif wall == 2
        return x1 * ymin - xmax * y1
    elseif wall == 3
        return x1 * ymin - xmin * y1
    else
        return x1 * ymax - xmin * y1
    end
end

"""
    connect_edges(x1, y1, x2, y2, wall1, wall2, xmin, xmax, ymin, ymax)

Calculates needed area component for shoelace formula between point one and
point two along the grid cell walls inbetween them.
Inputs:
    x1      <Number> x-coordinate of first point on wall 1
    y1      <Number> y-coordinate of first point on wall 1
    x2      <Number> x-coordinate of second point on wall 2
    y2      <Number> y-coordinate of second point on wall 2
    wall1   <Int> wall index (1-4) of point 1
    wall2   <Int> wall index (1-4) of point 2
    xmin    <AbstractFloat> minimum x-coordinate of grid cell
    xmax    <AbstractFloat> maximum x-coordinate of grid cell
    ymin    <AbstractFloat> minimum y-coordinate of grid cell
    ymax    <AbstractFloat> maximum y-coordinate of grid cell
Outputs:
    Area component of shoelace formula coming from the distance between point 1
    and point 2 along grid cell walls
"""
function connect_edges(x1::FT, y1, x2, y2, wall1, wall2,
    xmin, xmax, ymin, ymax,
) where {FT}
    connect_area = FT(0)
    # distance from point 1 to point two on same wall
    if wall1 == wall2 && is_clockwise_from(x1, y1, x2, y2, wall1)
        connect_area += x1 * y2 - x2 * y1
    else  # points are on a different wall, or point two is "behind" point one
        # From the point to the corner of wall 1
        connect_area += partial_edge_out(x1, y1, xmin, xmax, ymin, ymax, wall1)
        # Any intermediate walls (full length)
        if wall2 > wall1  # doesn't wrap around back past the north wall (1)
            for i in (wall1 + 1):(wall2 - 1)
                connect_area += full_edge_area(xmin, xmax, ymin, ymax, i)
            end
        else  # wraps back around the walls to or past the north wall (1)
            for i in (wall1 + 1):4
                connect_area += full_edge_area(xmin, xmax, ymin, ymax, i)
            end
            for i in 1:(wall2 - 1)
                connect_area += full_edge_area(xmin, xmax, ymin, ymax, i)
            end
        end
        # From the corner of wall 2 to the point
        connect_area += partial_edge_in(x2, y2, xmin, xmax, ymin, ymax, wall2)
    end
    return connect_area
end
