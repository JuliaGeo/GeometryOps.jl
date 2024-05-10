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
    intersection(geom_a, geom_b, [T::Type]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the intersection between two geometries as a list of geometries. Return an empty list
if none are found. The type of the list will be constrained as much as possible given the
input geometries. Furthermore, the user can provide a `target` type as a keyword argument and
a list of target geometries found in the intersection will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. If the
user is taking a intersection involving one or more multipolygons, and the multipolygon
might be comprised of polygons that intersect, if `fix_multipoly` is set to an
`IntersectingPolygons` correction (the default is `UnionIntersectingPolygons()`), then the
needed multipolygons will be fixed to be valid before performing the intersection to ensure
a correct answer. Only set `fix_multipoly` to nothing if you know that the multipolygons are
valid, as it will avoid unneeded computation. 

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
inter_points = GO.intersection(line1, line2; target = GI.PointTrait())
GI.coordinates.(inter_points)

# output
1-element Vector{Vector{Float64}}:
 [125.58375366067547, -14.83572303404496]
```
"""
function intersection(
    geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...,
) where {T<:AbstractFloat}
    return _intersection(TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b; kwargs...)
end

# Curve-Curve Intersections with target Point
_intersection(
    ::TraitTarget{GI.PointTrait}, ::Type{T},
    trait_a::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_a,
    trait_b::Union{GI.LineTrait, GI.LineStringTrait, GI.LinearRingTrait}, geom_b;
    kwargs...,
) where T = _intersection_points(T, trait_a, geom_a, trait_b, geom_b)

#= Polygon-Polygon Intersections with target Polygon
The algorithm to determine the intersection was adapted from "Efficient clipping
of efficient polygons," by Greiner and Hormann (1998).
DOI: https://doi.org/10.1145/274363.274364 =#
function _intersection(
    ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where {T}
    # First we get the exteriors of 'poly_a' and 'poly_b'
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then we find the intersection of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b, _inter_delay_cross_f, _inter_delay_bounce_f)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, _inter_step)
    if isempty(polys) # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b)
        if a_in_b
            push!(polys, GI.Polygon([tuples(ext_a)]))
        elseif b_in_a
            push!(polys, GI.Polygon([tuples(ext_b)]))
        end
    end
    remove_idx = falses(length(polys))
    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        hole_iterator = Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b)))
        _add_holes_to_polys!(T, polys, hole_iterator, remove_idx)
    end
    # Remove uneeded collinear points on same edge
    for p in polys
        _remove_collinear_points!(p, remove_idx)
    end
    return polys
end

# # Helper functions for Intersections with Greiner and Hormann Polygon Clipping

#= When marking the crossing status of a delayed crossing, the chain start point is bouncing
when the start point is a entry point and is a crossing point when the start point is an
exit point. The end of the chain has the opposite crossing / bouncing status. x is the 
entry/exit status. =#
_inter_delay_cross_f(x) = (!x, x)
#= When marking the crossing status of a delayed bouncing, the chain start and end points
are crossing if the current polygon's adjacent edges are within the non-tracing polygon. If
the edges are outside then the chain endpoints are marked as bouncing. x is a boolean
representing if the edges are inside or outside of the polygon. =#
_inter_delay_bounce_f(x, _) = x
#= When tracing polygons, step forward if the most recent intersection point was an entry
point, else step backwards where x is the entry/exit status. =#
_inter_step(x, _) =  x ? 1 : (-1)

#= Polygon with multipolygon intersection - note that all intersection regions between
`poly_a` and any of the sub-polygons of `multipoly_b` are counted as intersection polygons.
Unless specified with `fix_multipoly = nothing`, `multipolygon_b` will be validated using
the given (default is `UnionIntersectingPolygons()`) correction. =#
function _intersection(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent duplicated intersection regions
        multipoly_b = fix_multipoly(multipoly_b)
    end
    polys = Vector{_get_poly_type(T)}()
    for poly_b in GI.getpolygon(multipoly_b)
        append!(polys, intersection(poly_a, poly_b; target))
    end
    return polys
end

#= Multipolygon with polygon intersection is equivalent to taking the intersection of the
poylgon with the multipolygon and thus simply switches the order of operations and calls the
above method. =#
_intersection(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T = intersection(poly_b, multipoly_a; target , kwargs...)

#= Multipolygon with multipolygon intersection - note that all intersection regions between
any sub-polygons of `multipoly_a` and any of the sub-polygons of `multipoly_b` are counted
as intersection polygons. Unless specified with `fix_multipoly = nothing`, both 
`multipolygon_a` and `multipolygon_b` will be validated using the given (default is
`UnionIntersectingPolygons()`) correction. =#
function _intersection(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix both multipolygons to prevent duplicated regions
        multipoly_a = fix_multipoly(multipoly_a)
        multipoly_b = fix_multipoly(multipoly_b)
        fix_multipoly = nothing
    end
    polys = Vector{_get_poly_type(T)}()
    for poly_a in GI.getpolygon(multipoly_a)
        append!(polys, intersection(poly_a, multipoly_b; target, fix_multipoly))
    end
    return polys
end

# Many type and target combos aren't implemented
function _intersection(
    ::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b;
    kwargs...,
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

#= Calculates the intersection points between two lines if they exists and the fractional
component of each line from the initial end point to the intersection point where α is the
fraction along (a1, a2) and β is the fraction along (b1, b2).

Note that the first return is the type of intersection (line_cross, line_hinge, line_over,
or line_out). The type of intersection determines how many intersection points there are.
If the intersection is line_out, then there are no intersection points and the two
intersections aren't valid and shouldn't be used. If the intersection is line_cross or
line_hinge then the lines meet at one point and the first intersection is valid, while the
second isn't. Finally, if the intersection is line_over, then both points are valid and they
are the two points that define the endpoints of the overlapping region between the two
lines.

Also note again that each intersection is a tuple of two tuples. The first is the
intersection point (x,y) while the second is the ratio along the initial lines (α, β) for
that point. 

Calculation derivation can be found here: https://stackoverflow.com/questions/563198/ =#
function _intersection_point2(::Type{T}, (a1, a2)::Edge, (b1, b2)::Edge) where T
    # Default answer for no intersection
    line_orient = line_out
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    no_intr_result = (line_orient, intr1, intr2)
    # Seperate out line segment points
    (a1x, a1y), (a2x, a2y) = _tuple_point(a1, T), _tuple_point(a2, T)
    (b1x, b1y), (b2x, b2y) = _tuple_point(b1, T), _tuple_point(b2, T)
    # Check if envalopes of lines intersect
    a_ext = Extent(X = minmax(a1x, a2x), Y = minmax(a1y, a2y))
    b_ext = Extent(X = minmax(b1x, b2x), Y = minmax(b1y, b2y))
    !Extents.intersects(a_ext, b_ext) && return no_intr_result
    # Check orientation of two line segments with respect to one another
    a1_orient = Predicates.orient(b1, b2, a1)
    a2_orient = Predicates.orient(b1, b2, a2)
    a1_orient != 0 && a1_orient == a2_orient && return no_intr_result  # α < 0 or α > 1
    b1_orient = Predicates.orient(a1, a2, b1)
    b2_orient = Predicates.orient(a1, a2, b2)
    b1_orient != 0 && b1_orient == b2_orient && return no_intr_result  # β < 0 or β > 1
    # Intersection is collinear if all endpoints lie on the same line
    if a1_orient == a2_orient == b1_orient == b2_orient == 0
        line_orient, intr1, intr2 = _find_collinear_intersection(T, a1, a2, b1, b2, a_ext, b_ext)
    else  # Lines are not collinear and intersect in one point
        intr1 = if a1_orient == 0 || a2_orient == 0 || b1_orient == 0 || b2_orient == 0
            line_orient = line_hinge
            _find_hinge_intersection(T, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient)
        else
            line_orient = line_cross
            _find_cross_intersection(T, T, a1, a2, b1, b2, a_ext, b_ext)
        end
    end
    return line_orient, intr1, intr2
end

function _find_cross_intersection(::Type{T}, ::Type{P}, a1, a2, b1, b2, a_ext, b_ext) where {T, P}
    # First line runs from a to a + Δa
    (a1x, a1y), (a2x, a2y) = _tuple_point(a1, T), _tuple_point(a2, T)
    Δax, Δay = a2x - a1x, a2y - a1y
    # Second line runs from b to b + Δb 
    (b1x, b1y), (b2x, b2y) = _tuple_point(b1, T), _tuple_point(b2, T)
    Δbx, Δby = b2x - b1x, b2y - b1y
    # Differences between starting points
    Δba_x = b2x - b1x
    Δba_y = b2y - b1y
    Δa_cross_Δb = Δax * Δby - Δay * Δbx
    # Determine α value where 0 < α < 1
    α = T((Δba_x * Δby - Δba_y * Δbx) / Δa_cross_Δb)
    α = clamp(α, eps(T), one(T) - eps(T))
    # Determine β value where 0 < β < 1
    β = T((Δba_x * Δay - Δba_y * Δax) / Δa_cross_Δb)
    β = clamp(β, eps(T), one(T) - eps(T))
    # 
    x = T(a1x + α * Δax )#+ b1x + β * Δbx) / 2
    y = T(a1y + α * Δay )#+ b1y + β * Δby) / 2
    pt = (x, y)
    #
    # if equals(pt, a1)
    #     α_min = max(2eps(a1x) / Δax, 2eps(a1y) / Δay)
    #     @assert α_min > α "equals(pt, a1): $α $a1, $a2"
    #     x = a1x + α_min * Δax
    #     y = a1y + α_min * Δax
    #     α = α_min
    # elseif equals(pt, a2)
    #     α_min = max(2eps(a2x) / Δax, 2eps(a2y) / Δay)
    #     α_min = 1 - α_min
    #     @assert α_min < α "equals(pt, a2): $α $a1, $a2"
    #     x = a1x + α_min * Δax
    #     y = a1y + α_min * Δax
    #     α = α_min
    # elseif equals(pt, b1)
    #     β_min = max(2eps(b1x) / Δbx, 2eps(b1y) / Δby)
    #     @assert β_min > β "equals(pt, b1): $β $b1, $b2"
    #     x = b1x + β_min * Δbx
    #     y = b1y + β_min * Δbx
    #     β = β_min
    # elseif equals(pt, b2)
    #     β_min = max(2eps(b2x) / Δbx, 2eps(b2y) / Δby)
    #     β_min = 1 - β_min
    #     @assert β_min < β "equals(pt, b2): $β $b1, $b2"
    #     x = b1x + β_min * Δbx
    #     y = b1y + β_min * Δbx
    #     β = β_min
    # end
    return (x, y), (α, β)
end


# function _find_cross_intersection(::Type{T}, ::Type{P}, a1, a2, b1, b2, a_ext, b_ext) where {T, P}
#     (a1x, a1y), (a2x, a2y) = _tuple_point(a1, P), _tuple_point(a2, P)
#     (b1x, b1y), (b2x, b2y) = _tuple_point(b1, P), _tuple_point(b2, P)
#     a_dist, b_dist = distance(a1, a2, T), distance(b1, b2, T)
#     # Calculate needed values for intersection calculation
#     Δax = a1x - a2x
#     Δay = a1y - a2y
#     a_cross = (a1x * a2y) - (a2x * a1y)
#     Δbx = b1x - b2x
#     Δby = b1y - b2y
#     b_cross = (b1x * b2y) - (b2x * b1y)
#     det = (Δay * Δbx) - (Δby * Δax)
#     # Find x and y with given precision type P
#     x_p = ((Δax * b_cross) - (Δbx * a_cross)) / det
#     y_p = ((Δby * a_cross) - (Δay * b_cross)) / det
#     # Convert x and y to return type T
#     x, y = T(x_p), T(y_p)
#     pt = (x, y)
#     # Find α and β values for determined intersection point
#     α, β = distance(pt, a1, T) / a_dist, distance(pt, b1, T) / b_dist
#     # If point is outside line extents, adjust to nearest endpoint ± 2ϵ so point is in line
#     invalid_pt = !_point_in_extent(pt, a_ext) || !_point_in_extent(pt, b_ext) 
#     invalid_pt |= (equals(pt, a1) || equals(pt, a2) || equals(pt, b1) || equals(pt, b2))
#     invalid_pt |= (α ≤ 0 || α ≥ 1 || β ≤ 0 || β ≥ 1)
#     if invalid_pt
#         pt, α, β = _adjust_crossing_intersection(T, pt, α, β, a1, a2, b1, b2, a_ext, b_ext, a_dist, b_dist)
#     end
#     return pt, (α, β)
# end

#= Due to inexact floating point calculations, point can be outside of line extents. A
good approximation of the intersection is the endpoint of either a or b that is the
closest to the other segment. This might happen in extreem cases where the slopes of the
two segments are almost parallel, or where the endpoint of one segment lines almost on
another segment. However, since we know the lines cross, and not hinge, we will move
"in" on the segment by `2eps(end_point)` to make sure that we are within the line.

Another potential problem is that the intersection point is exactly equal to the
endpoint. In this case, again, the new intersection point should be that endpoint moved
"in" on the segment by `2eps(end_point)`.
=# 
function _adjust_crossing_intersection(::Type{T}, pt, α, β, a1, a2, b1, b2, a_ext, b_ext, a_dist, b_dist) where T
    a_line = GI.Line(StaticArrays.SVector(a1, a2))
    b_line = GI.Line(StaticArrays.SVector(b1, b2))
    local nearest_pt, nearest_α, nearest_β,  min_dist
    for (i, e) in enumerate((a1, a2, b1, b2))
        te = _tuple_point(e, T)
        ϵ_dist = sqrt(sum(2 .* eps.(te)).^2)
        ϵ_frac = ϵ_dist / (i < 3 ? a_dist : b_dist)
        ϵ_frac = max(ϵ_frac, 2eps(T))
        if iseven(i)
            ϵ_frac = 1 - ϵ_frac
            ϵ_frac = ϵ_frac == 1 ? prevfloat(one(T), 2) : ϵ_frac
        end
        near_pt, near_α, near_β, dist = if i < 3
            pt = a1 .+ (ϵ_frac .* (a2 .- a1))
            α = ϵ_frac
            β = distance(pt, b1, T) / b_dist
            d = distance(pt, b_line, T)
            pt, α, β, d
        else
            pt = b1 .+ (ϵ_frac .* (b2 .- b1))
            α = distance(pt, a1, T) / a_dist
            β = ϵ_frac
            d = distance(pt, a_line, T)
            pt, α, β, d
        end
        if i == 1 || dist < min_dist
            nearest_pt, nearest_α, nearest_β, min_dist = near_pt,  near_α, near_β, dist
        end
    end
    return nearest_pt, nearest_α, nearest_β
end

function _find_hinge_intersection(::Type{T}, a1, a2, b1, b2, a1_orient, a2_orient, b1_orient) where T
    pt, α, β = if equals(a1, b1)
        _tuple_point(a1, T), zero(T), zero(T)
    elseif equals(a1, b2)
        _tuple_point(a1, T), zero(T), one(T)
    elseif equals(a2, b1)
        _tuple_point(a2, T), one(T), zero(T)
    elseif equals(a2, b2)
        _tuple_point(a2, T), one(T), one(T)
    elseif a1_orient == 0
        β_val = _clamped_frac(distance(b1, a1, T), distance(b1, b2, T), eps(T))
        _tuple_point(a1, T), zero(T), β_val
    elseif a2_orient == 0
        β_val = _clamped_frac(distance(b1, a2, T), distance(b1, b2, T), eps(T))
        _tuple_point(a2, T), one(T), β_val
    elseif b1_orient == 0
        α_val = _clamped_frac(distance(a1, b1, T), distance(a1, a2, T), eps(T))
        _tuple_point(b1, T), α_val, zero(T)
    else  # b2_orient == 0
        α_val = _clamped_frac(distance(a1, b2, T), distance(a1, a2, T), eps(T))
        _tuple_point(b2, T), α_val, one(T)
    end
    return pt, (α, β)
end

function _find_collinear_intersection(::Type{T}, a1, a2, b1, b2, a_ext, b_ext) where T
    # Define default return for no intersection points
    line_orient = line_out
    intr1 = (zero(T), zero(T)), zero(T), zero(T)
    intr2 = intr1
    # Determine collinear line overlaps
    a1_in_b = _point_in_extent(a1, b_ext)
    a2_in_b = _point_in_extent(a2, b_ext)
    b1_in_a = _point_in_extent(b1, a_ext)
    b2_in_a = _point_in_extent(b2, a_ext)
    # Determine line distances
    a_dist, b_dist = distance(a1, a2, T), distance(b1, b2, T)
    _set_ab_collinear_intrs(args...) = _set_ab_collinear_intrs(args..., a1, b1, a_dist, b_dist)
    # Set collinear intersection points if they exist
    if a1_in_b && a2_in_b      # 1st vertex of a and 2nd vertex of a form overlap
        line_orient = line_over
        β1 = _clamped_frac(distance(a1, b1, T), b_dist)
        β2 = _clamped_frac(distance(a2, b1, T), b_dist)
        intr1 = (_tuple_point(a1, T), (zero(T), β1))
        intr2 = (_tuple_point(a2, T), (one(T), β2))
    elseif b1_in_a && b2_in_a  # 1st vertex of b and 2nd vertex of b form overlap
        line_orient = line_over
        α1 = _clamped_frac(distance(b1, a1, T), a_dist)
        α2 = _clamped_frac(distance(b2, a1, T), a_dist)
        intr1 = (_tuple_point(b1, T), (α1, zero(T)))
        intr2 = (_tuple_point(b2, T), (α2, one(T)))
    elseif a1_in_b && b1_in_a  # 1st vertex of a and 1st vertex of b form overlap
        if equals(a1, b1)
            line_orient = line_hinge
            intr1 = (_tuple_point(a1, T), zero(T), zero(T))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(a1, b1, zero(T), zero(T))
        end
    elseif a1_in_b && b2_in_a  # 1st vertex of a and 2nd vertex of b form overlap
        if equals(a1, b2)
            line_orient = line_hinge
            intr1 = (_tuple_point(a1, T), zero(T), one(T))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(a1, b2, zero(T), one(T)) 
        end
    elseif a2_in_b && b1_in_a  # 2nd vertex of a and 1st vertex of b form overlap
        if equals(a2, b1)
            line_orient = line_hinge
            intr1 = (_tuple_point(a2, T), one(T), zero(T))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(a2, b1, one(T), zero(T))
        end
    elseif a2_in_b && b2_in_a  # 2nd vertex of a and 2nd vertex of b form overlap
        if equals(a2, b2)
            line_orient = line_hinge
            intr1 = (_tuple_point(a2, T), one(T), one(T))
        else
            line_orient = line_over
            intr1, intr2 = _set_ab_collinear_intrs(a2, b2, one(T), one(T))
        end
    end
    return line_orient, intr1, intr2
end

_clamped_frac(x::T, y::T, ϵ = zero(T)) where T = clamp(x / y, ϵ, one(T) - ϵ)

_set_ab_collinear_intrs(::Type{T}, a_pt, b_pt, a_pt_α, b_pt_β, a1, b1, a_dist, b_dist) where T =
    (
        (_tuple_point(a_pt, T), (a_pt_α, _clamped_frac(distance(a_pt, b1, T), b_dist))),
        (_tuple_point(b_pt, T), (_clamped_frac(distance(b_pt, a1, T), a_dist), b_pt_β))
    )


function _intersection_point(::Type{T}, (a1, a2)::Edge, (b1, b2)::Edge) where T
    # Return line orientation and 2 intersection points + fractions (nothing if don't exist)
    line_orient = line_out
    intr1 = ((zero(T), zero(T)), (zero(T), zero(T)))
    intr2 = intr1
    # TODO: add an envelope check 
    # First line runs from p to p + r
    px, py = GI.x(a1), GI.y(a1)
    rx, ry = GI.x(a2) - px, GI.y(a2) - py
    # Second line runs from q to q + s 
    qx, qy = GI.x(b1), GI.y(b1)
    sx, sy = GI.x(b2) - qx, GI.y(b2) - qy

    # @show px + (rx * (((qx - px) * sy - (qy - py) * sx)/(rx * sy - ry * sx)))
    # @show (-px*ry*sx + rx*qx*sy - rx*qy*sx + rx*py*sx) / (rx * sy - ry * sx)
    # @show qx + (sx * (((qx - px) * ry - (qy - py) * rx)/(rx * sy - ry * sx)))
    # @show (qx*rx*sy - sx*px*ry - sx*qy*rx + sx*py*rx) / (rx * sy - ry * sx)
    # @show (qy - py + (ry/rx) * px - (sy/sx) * qx) / (ry/rx - sy/sx)
    # @show (rx*sx*(qy - py) + (ry * sx * px) - (sy * rx * qx)) / (ry * sx - sy * rx)
    # Intersections will be where p + αr = q + βs where 0 < α, β < 1
    Δqp_x = qx - px
    Δqp_y = qy - py
    if Predicates.isparallel((rx, ry), (sx, sy)) != 0  # non-parallel lines
        # Calculate α ratio if lines cross or touch
        a1_orient = Predicates.orient(b1, b2, a1)
        a2_orient = Predicates.orient(b1, b2, a2)
        # Lines don't cross α < 0 or α > 1
        a1_orient != 0 && a1_orient == a2_orient && return (line_orient, intr1, intr2)
        # Determine α value
        r_cross_s = rx * sy - ry * sx
        α, pt = if a1_orient == 0  # α = 0
            zero(T), (T(px), T(py))
        elseif a2_orient == 0  # α = 1
            one(T), (T(GI.x(a2)), T(GI.y(a2)))
        else # 0 < α < 1
            α_val = T((Δqp_x * sy - Δqp_y * sx) / r_cross_s)
            α_val = clamp(α_val, zero(T), one(T))
            α_val, (T(px + α_val * rx),  T(py + α_val * ry))
        end
        # Calculate β ratio if lines touch or cross
        b1_orient = Predicates.orient(a1, a2, b1)
        b2_orient = Predicates.orient(a1, a2, b2)
        # Lines don't cross β < 0 or β > 1
        b1_orient != 0 && b1_orient == b2_orient && return (line_orient, intr1, intr2)
        β, pt = if b1_orient == 0  # β = 0
            zero(T), (T(qx), T(qy))
        elseif b2_orient == 0  # β = 1
            one(T), (T(GI.x(b2)), T(GI.y(b2)))
        else  # 0 < β < 1
            β_val = T((Δqp_x * ry - Δqp_y * rx) / r_cross_s)
            β_val = clamp(β_val, zero(T), one(T))
            #= Floating point limitations could make intersection be endpoint if α≈0 or α≈1.
            In this case, see if multiplication by β gives a distinct number. Otherwise,
            replace with closest floating point number to endpoint.=#
            if (α != 0 && equals(a1, pt)) || (α != 1 && equals(a2, pt)) || equals(b1, pt) || equals(b2, pt)
                pt = (T(qx + β_val * sx), T(qy + β_val * sy))
                if equals(a1, pt)
                    α_min = max(eps(px) / rx, eps(py) / ry)
                    pt = (T(px + α_min * rx), T(py + α_min * ry))
                elseif equals(a2, pt)
                    α_max = 1 - max(eps(GI.x(a2)) / rx, eps(GI.y(a2)) / ry)
                    pt = (T(px + α_max * rx), T(py + α_max * ry))
                elseif equals(b1, pt)
                    β_min = max(eps(qx) / sx, eps(qy) / sy)
                    pt = (T(qx + β_min * sx), T(qy + β_min * sy))
                elseif equals(b2, pt)
                    β_max = 1 - max(eps(GI.x(b2)) / sx, eps(GI.y(b2)) / sy)
                    pt = (T(qx + β_max * sx), T(qy + β_max * sy))
                end
            end
            β_val, pt
        end
        # @show (T(qx + β * sx), T(qy + β * sy))
        # @show (T(px + α * rx),  T(py + α * ry))
        # Calculate intersection point using α and β
        # x, y = T(px + α * rx),  T(py + α * ry)
        # if (x == px && y == py) || (x == GI.x(a2) && y == GI.y(a2))
        #     x, y =  T(qx + β * sx),  T(qy + β * sy)
        # end
        intr1 = pt, (α, β)
        line_orient = (α == 0 || α == 1 || β == 0 || β == 1) ? line_hinge : line_cross
    elseif Predicates.iscollinear((Δqp_x, Δqp_y), (sx, sy)) == 0 # collinear parallel lines
        # Determine if lines touch or overlap and with what α and β values
        a1_side = Predicates.sameside(a1, b1, b2)
        a2_side = Predicates.sameside(a2, b1, b2)
        b1_side = Predicates.sameside(b1, a1, a2)
        b2_side = Predicates.sameside(b2, a1, a2)
        # Lines touch or overlap if endpoints of line a are on/in line b and visa versa
        r_dot_s = rx * sx + ry * sy
        # Determine which endpoints start and end the overlapping region
        n_intrs = 0
        if a1_side != 1 || a2_side != 1  # at least one endpoint of line a is in/on line b
            s_dot_s = sx^2 + sy^2
            a1_β = T(-(Δqp_x * sx + Δqp_y * sy) / s_dot_s)
            if a1_side != 1  # 0 ≤ a1_β ≤ 1
                n_intrs += 1
                a1_β = if a1_side == 0  # a1_β == 0 or  a1_β == 1
                    equals(a1, b1) ? zero(T) : one(T)
                else  # 0 < a1_β < 1
                    clamp(a1_β, zero(T), one(T))
                end
                intr1 = (T.(a1), (zero(T), a1_β))
            end
            if a2_side != 1  # 0 ≤ a2_β ≤ 1
                n_intrs += 1
                a2_β = if a2_side == 0  # a2_β == 0 or  a2_β == 1
                    equals(a2, b1) ? zero(T) : one(T)
                else  # 0 < a2_β < 1
                    β_val = a1_β + r_dot_s / s_dot_s
                    clamp(T(β_val), zero(T), one(T))
                end
                new_intr = (T.(a2), (one(T), a2_β))
                n_intrs == 1 && (intr1 = new_intr)
                n_intrs == 2 && (intr2 = new_intr)
            end
        end
        if b1_side == -1 || b2_side == -1  # at least one endpoint of line b is in line a
            r_dot_r = (rx^2 + ry^2)
            b1_α = T((Δqp_x * rx + Δqp_y * ry) / r_dot_r)
            if b1_side == -1   # 0 < b1_α < 1
                n_intrs += 1
                b1_α = clamp(b1_α, zero(T), one(T))
                new_intr = (T.(b1), (b1_α, zero(T)))
                n_intrs == 1 && (intr1 = new_intr)
                n_intrs == 2 && (intr2 = new_intr)
            end
            if b2_side == -1  # 0 < b2_α < 1
                n_intrs += 1
                b2_α = T(b1_α + r_dot_s / r_dot_r)
                b2_α = clamp(b2_α, zero(T), one(T))
                new_intr = (T.(b2), (b2_α, one(T)))
                n_intrs == 1 && (intr1 = new_intr)
                n_intrs == 2 && (intr2 = new_intr)
            end
        end
        if n_intrs == 1
            line_orient = line_hinge
        elseif n_intrs > 1
            line_orient = line_over
        end
    end
    return line_orient, intr1, intr2
end