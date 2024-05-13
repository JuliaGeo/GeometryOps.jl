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
    return _intersection(
        TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b;
        exact = _True(), kwargs...,
    )
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
    exact, kwargs...,
) where {T}
    # First we get the exteriors of 'poly_a' and 'poly_b'
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then we find the intersection of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b, _inter_delay_cross_f, _inter_delay_bounce_f; exact)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, _inter_step)
    if isempty(polys) # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b; exact)
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
        _add_holes_to_polys!(T, polys, hole_iterator, remove_idx; exact)
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
function _intersection_points(::Type{T}, ::GI.AbstractTrait, a, ::GI.AbstractTrait, b; exact = _False()) where T
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
            line_orient, intr1, _ = _intersection_point(T, edges_a[i], edges_b[j]; exact)
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
function _intersection_point(::Type{T}, (a1, a2)::Edge, (b1, b2)::Edge; exact) where T
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
    # Intersections will be where p + αr = q + βs where 0 < α, β < 1
    Δqp_x = qx - px
    Δqp_y = qy - py
    r_cross_s = rx * sy - ry * sx
    if Predicates.cross((rx, ry), (sx, sy); exact) != 0  # non-parallel lines
        # Calculate α ratio if lines cross or touch
        a1_orient = Predicates.orient(b1, b2, a1)
        a2_orient = Predicates.orient(b1, b2, a2)
        # Lines don't cross α < 0 or α > 1
        a1_orient != 0 && a1_orient == a2_orient && return (line_orient, intr1, intr2)
        # Determine α value
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
        # Calculate intersection point using α and β
        # x, y = T(px + α * rx),  T(py + α * ry)
        # if (x == px && y == py) || (x == GI.x(a2) && y == GI.y(a2))
        #     x, y =  T(qx + β * sx),  T(qy + β * sy)
        # end
        intr1 = pt, (α, β)
        line_orient = (α == 0 || α == 1 || β == 0 || β == 1) ? line_hinge : line_cross
    elseif Predicates.cross((Δqp_x, Δqp_y), (sx, sy); exact) == 0 # collinear parallel lines
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