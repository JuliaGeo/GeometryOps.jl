# # Difference Polygon Clipping
export difference


"""
    difference(geom_a, geom_b, [T::Type]; target::Type)

Return the difference between two geometries as a list of geometries. Return an empty list
if none are found. The type of the list will be constrained as much as possible given the
input geometries. Furthermore, the user can provide a `taget` type as a keyword argument and
a list of target geometries found in the difference will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. 

## Example 

```jldoctest
import GeoInterface as GI, GeometryOps as GO

poly1 = GI.Polygon([[[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]])
poly2 = GI.Polygon([[[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]]])
diff_poly = GO.difference(poly1, poly2; target = GI.PolygonTrait())
GI.coordinates.(diff_poly)

# output
1-element Vector{Vector{Vector{Vector{Float64}}}}:
 [[[6.5, 3.5], [5.0, 5.0], [0.0, 0.0], [5.0, -5.0], [6.5, -3.5], [3.0, 0.0], [6.5, 3.5]]]
```
"""
function difference(
    geom_a, geom_b, ::Type{T} = Float64; target=nothing, kwargs...,
) where {T<:AbstractFloat}
    return _difference(TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b; kwargs...)
end

#= The 'difference' function returns the difference of two polygons as a list of polygons.
The algorithm to determine the difference was adapted from "Efficient clipping of efficient
polygons," by Greiner and Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _difference(
    ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...
) where T
    # Get the exterior of the polygons
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Find the difference of the exterior of the polygons
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b, _diff_delay_cross_f, _diff_delay_bounce_f)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, _diff_step)
    # if no crossing points, determine if either poly is inside of the other
    if isempty(polys)
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b)
        # add case for if they polygons are the same (all intersection points!)
        # add a find_first check to find first non-inter poly!
        if b_in_a && !a_in_b  # b in a and can't be the same polygon
            poly_a_b_hole = GI.Polygon([tuples(ext_a), tuples(ext_b)])
            push!(polys, poly_a_b_hole)
        elseif !b_in_a && !a_in_b # polygons don't intersect
            push!(polys, tuples(poly_a))
            return polys
        end
    end
    remove_idx = falses(length(polys))
    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0
        _add_holes_to_polys!(T, polys, GI.gethole(poly_a), remove_idx)
    end
    if GI.nhole(poly_b) != 0
        for hole in GI.gethole(poly_b)
            hole_poly = GI.Polygon(StaticArrays.SVector(hole))
            new_polys = intersection(hole_poly, poly_a, T; target = GI.PolygonTrait)
            if length(new_polys) > 0
                append!(polys, new_polys)
            end
        end
    end
    # Remove uneeded collinear points on same edge
    for p in polys
        _remove_collinear_points!(p, remove_idx)
    end
    return polys
end

# # Helper functions for Differences with Greiner and Hormann Polygon Clipping

#= When marking the crossing status of a delayed crossing, the chain start point is crossing
when the start point is a entry point and is a bouncing point when the start point is an
exit point. The end of the chain has the opposite crossing / bouncing status. =#
_diff_delay_cross_f(x) = (x, !x)
#= When marking the crossing status of a delayed bouncing, the chain start and end points
are crossing if the current polygon's adjacent edges are within the non-tracing polygon and
we are tracing b_list or if the edges are outside and we are on a_list. Otherwise the
endpoints are marked as crossing. x is a boolean representing if the edges are inside or
outside of the polygon and y is a variable that is true if we are on a_list and false if we
are on b_list. =#
_diff_delay_bounce_f(x, y) = x ⊻ y
#= When tracing polygons, step forwards if the most recent intersection point was an entry
point and we are currently tracing b_list or if it was an exit point and we are currently
tracing a_list, else step backwards, where x is the entry/exit status and y is a variable
that is true if we are on a_list and false if we are on b_list. =#
_diff_step(x, y) = (x ⊻ y) ? 1 : (-1)

function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    kwargs...,
) where T
    polys = [tuples(poly_a, T)]
    for poly_b in GI.getpolygon(multipoly_b)
        isempty(polys) && break
        polys = mapreduce(p -> difference(p, poly_b; target = target), append!, polys)
    end
    return polys
end

function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    fix_multipoly = true, kwargs...,
) where T
    if fix_multipoly
        multipoly_a = MinimalMultiPolygon()(multipoly_a)
    end
    # TODO: Should we fix multipoly_a -> shouldn't happen every time if called with 2 multipoly
    polys = Vector{_get_poly_type(T)}()
    sizehint!(polys, GI.npolygon(multipoly_a))
    for poly_a in GI.getpolygon(multipoly_a)
        append!(polys, difference(poly_a, poly_b; target = target))
    end
    return polys
end

function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = true, kwargs...,
) where T
    if fix_multipoly
        multipoly_a = MinimalMultiPolygon()(multipoly_a)
        fix_multipoly = false
    end
    local polys
    for (i, poly_b) in enumerate(GI.getpolygon(multipoly_b))
        polys = difference(i == 1 ? multipoly_a : GI.MultiPolygon(polys), poly_b;
            target = target, fix_multipoly = fix_multipoly)
    end
    return polys
end

# Many type and target combos aren't implemented
function _difference(
    ::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
) where {Target, T}
    @assert(
        false,
        "Difference between $trait_a and $trait_b with target $Target isn't implemented yet.",
    )
    return nothing
end

