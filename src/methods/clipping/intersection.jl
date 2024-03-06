# # Geometry Intersection
export intersection, intersection_points

"""
    Enum LineOrientation
Enum for the orientation of a line with respect to a curve. A line can be
`line_cross` (crossing over the curve), `line_hinge` (crossing the endpoint of the curve),
`line_over` (colinear with the curve), or `line_out` (not interacting with the curve).
"""
@enum LineOrientation line_cross=1 line_hinge=2 line_over=3 line_out=4

"""
    intersection(geom_a, geom_b, [T::Type]; target::Type)

Return the intersection between two geometries as a list of geometries. Return an empty list
if none are found. The type of the list will be constrained as much as possible given the
input geometries. Furthermore, the user can provide a `taget` type as a keyword argument and
a list of target geometries found in the intersection will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. 

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
inter_points = GO.intersection(line1, line2; target = GI.PointTrait)
GI.coordinates.(inter_points)

# output
1-element Vector{Vector{Float64}}:
 [125.58375366067547, -14.83572303404496]
```
"""
function intersection(
    geom_a, geom_b, ::Type{T} = Float64; target::Type{Target} = Nothing,
) where {T <: AbstractFloat, Target <: Union{Nothing, GI.AbstractTrait}}
    return _intersection(Target, T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end

# Curve-Curve Intersections with target Point
_intersection(
    ::Type{GI.PointTrait}, ::Type{T},
    trait_a::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_a,
    trait_b::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_b,
) where T = _intersection_points(T, trait_a, geom_a, trait_b, geom_b)


#= Polygon-Polygon Intersections with target Polygon
The algorithm to determine the intersection was adapted from "Efficient clipping
of efficient polygons," by Greiner and Hormann (1998).
DOI: https://doi.org/10.1145/274363.274364 =#
function _intersection(
    ::Type{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b,
) where {T}
    # First we get the exteriors of 'poly_a' and 'poly_b'
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then we find the intersection of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, (x, y) -> x ? 1 : (-1))
    if isempty(polys) # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b)
        if a_in_b
            push!(polys, GI.Polygon([tuples(ext_a)]))
        elseif b_in_a
            push!(polys, GI.Polygon([tuples(ext_b)]))
        end
    end
    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        hole_iterator = Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b)))
        _add_holes_to_polys!(T, polys, hole_iterator)
    end    
    return polys
end

# Many type and target combos aren't implemented
function _intersection(
    ::Type{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
) where {Target, T}
    @assert(
        false,
        "Intersection between $trait_a and $trait_b with target $Target isn't implemented yet.",
    )
    return nothing
end

"""
    intersection_points(
        geom_a,
        geom_b,
    )::Union{
        ::Vector{::Tuple{::Real, ::Real}},
        ::Nothing,
    }

Return a list of intersection points between two geometries of type GI.Point.
If no intersection point was possible given geometry extents, returns an empty
list.
"""
intersection_points(geom_a, geom_b, ::Type{T} = Float64) where T <: AbstractFloat =
    _intersection_points(T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)


#= Calculates the list of intersection points between two geometries, inlcuding line
segments, line strings, linear rings, polygons, and multipolygons. If no intersection points
were possible given geometry extents or if none are found, return an empty list of
GI.Points. =#
function _intersection_points(::Type{T}, ::GI.AbstractTrait, a, ::GI.AbstractTrait, b) where T
    # Initialize an empty list of points
    result = GI.Point[]
    # Check if the geometries extents even overlap
    Extents.intersects(GI.extent(a), GI.extent(b)) || return result
    # Create a list of edges from the two input geometries
    edges_a, edges_b = map(sort! ∘ to_edges, (a, b))
    npoints_a, npoints_b  = length(edges_a), length(edges_b)
    a_closed = npoints_a > 1 && edges_a[1][1] == edges_a[end][1]
    b_closed = npoints_b > 1 && edges_b[1][1] == edges_b[end][1]
    if npoints_a > 0 && npoints_b > 0
        # Loop over pairs of edges and add any intersection points to results
        for i in eachindex(edges_a), j in eachindex(edges_b)
            line_orient, intr1, _ = _intersection_point(T, edges_a[i], edges_b[j])
            # TODO: Add in degenerate intersection points when line_over
            if line_orient == line_cross || line_orient == line_hinge
                #=
                Determine if point is on edge (all edge endpoints excluded
                except for the last edge for an open geometry)
                =#
                point, (α, β) = intr1
                on_a_edge = (!a_closed && i == npoints_a && 0 <= α <= 1) ||
                    (0 <= α < 1)
                on_b_edge = (!b_closed && j == npoints_b && 0 <= β <= 1) ||
                    (0 <= β < 1)
                if on_a_edge && on_b_edge
                    push!(result, GI.Point(point))
                end
            end
        end
    end
    return result
end

#= Calculates the intersection point between two lines if it exists, and as if the line
extended to infinity, and the fractional component of each line from the initial end point
to the intersection point.
Inputs:
    (a1, a2)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} first line
    (b1, b2)::Tuple{Tuple{::Real, ::Real}, Tuple{::Real, ::Real}} second line
Outputs:
    (x, y)::Tuple{::Real, ::Real} intersection point
    (t, u)::Tuple{::Real, ::Real} fractional length of lines to intersection
    Both are ::Nothing if point doesn't exist!

Calculation derivation can be found here:
    https://stackoverflow.com/questions/563198/
=#
function _intersection_point(::Type{T}, (a1, a2)::Tuple, (b1, b2)::Tuple) where T
    # Return line orientation and 2 intersection points + fractions (nothing if don't exist)
    line_orient = line_out
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    # First line runs from p to p + r
    px, py = GI.x(a1), GI.y(a1)
    rx, ry = GI.x(a2) - px, GI.y(a2) - py
    # Second line runs from q to q + s 
    qx, qy = GI.x(b1), GI.y(b1)
    sx, sy = GI.x(b2) - qx, GI.y(b2) - qy
    # Intersections will be where p + αr = q + βs where 0 < α, β < 1 and
    r_cross_s = rx * sy - ry * sx
    Δqp_x = qx - px
    Δqp_y = qy - py
    if r_cross_s != 0  # if lines aren't parallel
        α = (Δqp_x * sy - Δqp_y * sx) / r_cross_s
        β = (Δqp_x * ry - Δqp_y * rx) / r_cross_s
        x = px + α * rx
        y = py + α * ry
        if 0 ≤ α ≤ 1 && 0 ≤ β ≤ 1
            intr1 = (T(x), T(y)), (T(α), T(β))
            line_orient = (α == 0 || α == 1 || β == 0 || β == 1) ? line_hinge : line_cross
        end
    elseif sx * Δqp_y == sy * Δqp_x  # if parallel lines are collinear
        # Determine overlap fractions
        r_dot_r = (rx^2 + ry^2)
        s_dot_s = (sx^2 + sy^2)
        r_dot_s = rx * sx + ry * sy
        b1_α = (Δqp_x * rx + Δqp_y * ry) / r_dot_r
        b2_α = b1_α + r_dot_s / r_dot_r
        a1_β = -(Δqp_x * sx + Δqp_y * sy) / s_dot_s
        a2_β = a1_β + r_dot_s / s_dot_s
        # Determine which endpoints start and end the overlapping region
        n_intrs = 0
        if 0 ≤ a1_β ≤ 1
            n_intrs += 1
            intr1 = (T.(a1), (zero(T), T(a1_β)))
        end
        if 0 ≤ a2_β ≤ 1
            n_intrs += 1
            new_intr = (T.(a2), (one(T), T(a2_β)))
            n_intrs == 1 && (intr1 = new_intr)
            n_intrs == 2 && (intr2 = new_intr)
        end
        if 0 < b1_α < 1 
            n_intrs += 1
            new_intr = (T.(b1), (T(b1_α), zero(T)))
            n_intrs == 1 && (intr1 = new_intr)
            n_intrs == 2 && (intr2 = new_intr)
        end
        if 0 < b2_α < 1
            n_intrs += 1
            new_intr = (T.(b2), (T(b2_α), one(T)))
            n_intrs == 1 && (intr1 = new_intr)
            n_intrs == 2 && (intr2 = new_intr)
        end
        if n_intrs == 1
            line_orient = line_hinge
        elseif n_intrs == 2
            line_orient = line_over
        end
    end
    return line_orient, intr1, intr2
end
