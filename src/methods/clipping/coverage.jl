export coverage

#=
## What is coverage?

Coverage is the amount of geometry area within a bounding box defined by the minimum and
maximum x and y-coordinates of that bounding box, or an Extent containing that information.

To provide an example, consider this rectangle:
```@example rect
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

rect = GI.Polygon([[(-1,0), (-1,1), (1,1), (1,0), (-1,0)]])
cell = GI.Polygon([[(0, 0), (0, 2), (2, 2), (2, 0), (0, 0)]])
xmin, xmax, ymin, ymax = 0, 2, 0, 2
f, a, p = poly(collect(GI.getpoint(cell)); axis = (; aspect = DataAspect()))
poly!(collect(GI.getpoint(rect)))
f
```
It is clear that half of the polygon is within the cell, so the coverage should be 1.0, half
of the area of the rectangle. 
```@example rect
GO.coverage(rect, xmin, xmax, ymin, ymax)
```

## Implementation

This is the GeoInterface-compatible implementation. First, we implement a wrapper method
that dispatches to the correct implementation based on the geometry trait. This is also used
in the implementation, since it's a lot less work!

Note that the coverage is zero for all points and curves, even if the curves are closed like
with a linear ring.
=#

# Targets for applys functions
const _COVERAGE_TARGETS = TraitTarget{Union{GI.PolygonTrait,GI.AbstractCurveTrait,GI.MultiPointTrait,GI.PointTrait}}()

# Wall types for coverage - used to identify which cell boundary we're on
const UNKNOWN, NORTH, EAST, SOUTH, WEST = 0:4

"""
    coverage(geom, xmin, xmax, ymin, ymax, [T = Float64])::T

Returns the area of intersection between given geometry and grid cell defined by its minimum
and maximum x and y-values. This is computed differently for different geometries:

- The signed area of a point is always zero.
- The signed area of a curve is always zero.
- The signed area of a polygon is calculated by tracing along its edges and switching to the
    cell edges if needed.
- The coverage of a geometry collection, multi-geometry, feature collection of
    array/iterable is the sum of the coverages of all of the sub-geometries. 

Result will be of type T, where T is an optional argument with a default value
of Float64.
"""
function coverage(geom, xmin, xmax, ymin, ymax, ::Type{T} = Float64; threaded=false) where T <: AbstractFloat
    applyreduce(+, _COVERAGE_TARGETS, geom; threaded, init=zero(T)) do g
        _coverage(T, GI.trait(g), g, T(xmin), T(xmax), T(ymin), T(ymax))
    end
end

function coverage(geom, cell_ext::Extents.Extent, ::Type{T} = Float64; threaded=false) where T <: AbstractFloat
    (xmin, xmax), (ymin, ymax) = values(cell_ext)
    return coverage(geom, xmin, xmax, ymin, ymax, T; threaded = threaded)
end

# Points, MultiPoints, Curves, MultiCurves have zero coverage
_coverage(::Type{T}, ::GI.AbstractGeometryTrait, geom, xmin, xmax, ymin, ymax; kwargs...) where T = zero(T)

# Polygons: compute coverage for exterior and subtract coverage of holes
function _coverage(::Type{T}, ::GI.PolygonTrait, poly, xmin, xmax, ymin, ymax; exact = False()) where T
    GI.isempty(poly) && return zero(T)
    
    # Create a polygon representing the cell
    cell_poly = GI.Polygon([[(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]])
    
    # Compute intersection between input polygon and cell
    intersection_polys = intersection(poly, cell_poly; target = GI.PolygonTrait())
    
    # Sum up areas of all intersection polygons
    total_area = zero(T)
    for p in intersection_polys
        total_area += area(p)
    end
    
    return total_area
end

# Helper: Check if point is strictly inside cell (not on boundary)
function _point_strictly_inside(p, xmin, xmax, ymin, ymax)
    x, y = GI.x(p), GI.y(p)
    return xmin < x < xmax && ymin < y < ymax
end

# Helper: Check if point is on cell boundary and determine which wall
function _point_on_boundary(p, xmin, xmax, ymin, ymax)
    x, y = GI.x(p), GI.y(p)
    if x == xmin && ymin ≤ y ≤ ymax
        return WEST
    elseif x == xmax && ymin ≤ y ≤ ymax
        return EAST
    elseif y == ymin && xmin ≤ x ≤ xmax
        return SOUTH
    elseif y == ymax && xmin ≤ x ≤ xmax
        return NORTH
    else
        return UNKNOWN
    end
end

# Helper: Check if point is exactly at a vertex
function _point_at_vertex(p, xmin, xmax, ymin, ymax)
    x, y = GI.x(p), GI.y(p)
    return (x == xmin && (y == ymin || y == ymax)) ||
           (x == xmax && (y == ymin || y == ymax))
end

# Helper: Check if value b is between a and c (inclusive)
function _between(b, a, c)
    return min(a, c) ≤ b ≤ max(a, c)
end

# Helper: Find intersection of line segment with cell boundary
function _line_intersect_cell(::Type{T}, (x1, y1), (x2, y2), xmin, xmax, ymin, ymax) where T
    Δx, Δy = x2 - x1, y2 - y1
    
    # Handle vertical lines
    if Δx == 0
        if x1 < xmin || x1 > xmax
            return (UNKNOWN, (zero(T), zero(T))), (UNKNOWN, (zero(T), zero(T)))
        end
        t_min = (ymin - y1) / Δy
        t_max = (ymax - y1) / Δy
        if t_min > t_max
            t_min, t_max = t_max, t_min
        end
        if t_min > 1 || t_max < 0
            return (UNKNOWN, (zero(T), zero(T))), (UNKNOWN, (zero(T), zero(T)))
        end
        t_min = max(0, t_min)
        t_max = min(1, t_max)
        return (SOUTH, (x1, y1 + t_min * Δy)), (NORTH, (x1, y1 + t_max * Δy))
    end
    
    # Handle horizontal lines
    if Δy == 0
        if y1 < ymin || y1 > ymax
            return (UNKNOWN, (zero(T), zero(T))), (UNKNOWN, (zero(T), zero(T)))
        end
        t_min = (xmin - x1) / Δx
        t_max = (xmax - x1) / Δx
        if t_min > t_max
            t_min, t_max = t_max, t_min
        end
        if t_min > 1 || t_max < 0
            return (UNKNOWN, (zero(T), zero(T))), (UNKNOWN, (zero(T), zero(T)))
        end
        t_min = max(0, t_min)
        t_max = min(1, t_max)
        return (WEST, (x1 + t_min * Δx, y1)), (EAST, (x1 + t_max * Δx, y1))
    end
    
    # General case: parametric line equation
    tx1 = (xmin - x1) / Δx  # Intersection with x = xmin
    tx2 = (xmax - x1) / Δx  # Intersection with x = xmax
    ty1 = (ymin - y1) / Δy  # Intersection with y = ymin
    ty2 = (ymax - y1) / Δy  # Intersection with y = ymax
    
    # Find valid intersections (0 ≤ t ≤ 1)
    ts = Tuple{T,Int}[]
    for (t, wall) in [(tx1, WEST), (tx2, EAST), (ty1, SOUTH), (ty2, NORTH)]
        if 0 ≤ t ≤ 1
            x = x1 + t * Δx
            y = y1 + t * Δy
            if xmin ≤ x ≤ xmax && ymin ≤ y ≤ ymax
                push!(ts, (t, wall))
            end
        end
    end
    
    # Sort by parameter t to get intersections in order
    sort!(ts)
    
    if length(ts) < 2
        return (UNKNOWN, (zero(T), zero(T))), (UNKNOWN, (zero(T), zero(T)))
    end
    
    t1, wall1 = ts[1]
    t2, wall2 = ts[2]
    p1 = (x1 + t1 * Δx, y1 + t1 * Δy)
    p2 = (x1 + t2 * Δx, y1 + t2 * Δy)
    
    return (wall1, p1), (wall2, p2)
end

#= Calculates the area of the filled ring within the cell defined by corners with (xmin, ymin),
(xmin, ymax), (xmax, ymax), and (xmax, ymin). =#
function _coverage(::Type{T}, ring, xmin, xmax, ymin, ymax; exact) where T
    @info "Starting coverage calculation" xmin xmax ymin ymax
    cov_area = zero(T)
    
    # Get all vertices
    points = [_tuple_point(p, T) for p in GI.getpoint(ring)]
    n = length(points)
    n < 3 && return zero(T)  # Not a valid ring
    
    # Process each edge
    for i in 1:n
        p1 = points[i]
        p2 = points[i == n ? 1 : i + 1]
        @info "Processing edge" i p1 p2
        
        # Check point positions
        p1_inside = _point_strictly_inside(p1, xmin, xmax, ymin, ymax)
        p2_inside = _point_strictly_inside(p2, xmin, xmax, ymin, ymax)
        p1_wall = _point_on_boundary(p1, xmin, xmax, ymin, ymax)
        p2_wall = _point_on_boundary(p2, xmin, xmax, ymin, ymax)
        p1_vertex = _point_at_vertex(p1, xmin, xmax, ymin, ymax)
        p2_vertex = _point_at_vertex(p2, xmin, xmax, ymin, ymax)
        
        @info "Point classifications" p1_inside p2_inside p1_wall p2_wall p1_vertex p2_vertex
        
        # Skip if both points are the same vertex
        if p1_vertex && p2_vertex && p1 == p2
            @info "Skipping identical vertex points"
            continue
        end
        
        # Case 1: Both points inside - add full edge contribution
        if p1_inside && p2_inside
            area = _area_component(p1, p2)
            @info "Both points inside, adding area" area
            cov_area += area
            continue
        end
        
        # Case 2: Both points on boundary
        if p1_wall != UNKNOWN && p2_wall != UNKNOWN
            # Only add if not the same point
            if p1 != p2
                area = _area_component(p1, p2)
                @info "Both points on boundary, adding area" area p1_wall p2_wall
                cov_area += area
            else
                @info "Skipping identical boundary points"
            end
            continue
        end
        
        # Case 3: Edge crosses cell - find intersections
        inter1, inter2 = _line_intersect_cell(T, p1, p2, xmin, xmax, ymin, ymax)
        @info "Found intersections" inter1 inter2
        
        # No valid intersections
        if inter1[1] == UNKNOWN
            @info "No valid intersections"
            continue
        end
        
        # Determine the part of the edge that's inside the cell
        if p1_inside || p1_wall != UNKNOWN
            area = _area_component(p1, inter2[2])
            @info "P1 inside/on boundary, adding area" area
            cov_area += area
        elseif p2_inside || p2_wall != UNKNOWN
            area = _area_component(inter1[2], p2)
            @info "P2 inside/on boundary, adding area" area
            cov_area += area
        else
            area = _area_component(inter1[2], inter2[2])
            @info "Using intersection points, adding area" area
            cov_area += area
        end
    end
    
    # Return absolute area / 2 (shoelace formula)
    cov_area = abs(cov_area) / 2
    @info "Raw area" cov_area
    
    # Special case: if area is 0 but cell center is inside polygon
    if cov_area == 0
        center = ((xmin + xmax)/2, (ymin + ymax)/2)
        if _point_filled_curve_orientation(center, ring; in = true, on = false, out = false, exact)
            @info "Cell center inside polygon, using cell area"
            return (xmax - xmin) * (ymax - ymin)
        end
    end
    
    @info "Final area" cov_area
    return cov_area
end
