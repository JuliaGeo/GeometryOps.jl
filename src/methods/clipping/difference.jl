# # Difference Polygon Clipping
export difference


"""
    difference(geom_a, geom_b, [T::Type]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the difference between two geometries as a list of geometries. Return an empty list
if none are found. The type of the list will be constrained as much as possible given the
input geometries. Furthermore, the user can provide a `taget` type as a keyword argument and
a list of target geometries found in the difference will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. If the
user is taking a intersection involving one or more multipolygons, and the multipolygon
might be comprised of polygons that intersect, if `fix_multipoly` is set to an
`IntersectingPolygons` correction (the default is `UnionIntersectingPolygons()`), then the
needed multipolygons will be fixed to be valid before performing the intersection to ensure
a correct answer. Only set `fix_multipoly` to false if you know that the multipolygons are
valid, as it will avoid unneeded computation. 

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
    return _difference(
        TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b;
        exact = True(), kwargs...,
    )
end

#= The 'difference' function returns the difference of two polygons as a list of polygons.
The algorithm to determine the difference was adapted from "Efficient clipping of efficient
polygons," by Greiner and Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _difference(
    ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    exact, kwargs...
) where T
    # Get the exterior of the polygons
    ext_a = LazyClosedRing(GI.getexterior(poly_a))
    ext_b = LazyClosedRing(GI.getexterior(poly_b))
    # Find the difference of the exterior of the polygons
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b, _diff_delay_cross_f, _diff_delay_bounce_f; exact)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, _diff_step, poly_a, poly_b)
    # if no crossing points, determine if either poly is inside of the other
    if isempty(polys)
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b; exact)
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
        _add_holes_to_polys!(T, polys, GI.gethole(poly_a), remove_idx; exact)
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
    # Remove unneeded collinear points on same edge
    _remove_collinear_points!(polys, remove_idx, poly_a, poly_b)
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

#= Polygon with multipolygon difference - note that all intersection regions between
`poly_a` and any of the sub-polygons of `multipoly_b` are removed from `poly_a`. =#
function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    kwargs...,
) where T
    polys = [tuples(poly_a, T)]
    for poly_b in GI.getpolygon(multipoly_b)
        isempty(polys) && break
        polys = mapreduce(p -> difference(p, poly_b; target), append!, polys)
    end
    return polys
end

#= Multipolygon with polygon difference - note that all intersection regions between
sub-polygons of `multipoly_a` and `poly_b` will be removed from the corresponding
sub-polygon. Unless specified with `fix_multipoly = nothing`, `multipolygon_a` will be
validated using the given (default is `UnionIntersectingPolygons()`) correction. =#
function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_a to prevent returning an invalid multipolygon
        multipoly_a = fix_multipoly(multipoly_a)
    end
    polys = Vector{_get_poly_type(T)}()
    sizehint!(polys, GI.npolygon(multipoly_a))
    for poly_a in GI.getpolygon(multipoly_a)
        append!(polys, difference(poly_a, poly_b; target))
    end
    return polys
end

#= Multipolygon with multipolygon difference - note that all intersection regions between
sub-polygons of `multipoly_a` and sub-polygons of `multipoly_b` will be removed from the
corresponding sub-polygon of `multipoly_a`. Unless specified with `fix_multipoly = nothing`,
`multipolygon_a` will be validated using the given (default is `UnionIntersectingPolygons()`)
correction. =#
function _difference(
    target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_a to prevent returning an invalid multipolygon
        multipoly_a = fix_multipoly(multipoly_a)
        fix_multipoly = nothing
    end
    local polys
    for (i, poly_b) in enumerate(GI.getpolygon(multipoly_b))
        #= Removing intersections of `multipoly_a`` with pieces of `multipoly_b`` - as
        pieces of `multipolygon_a`` are removed, continue to take difference with new shape
        `polys` =#
        polys = if i == 1
            difference(multipoly_a, poly_b; target, fix_multipoly)
        else
            difference(GI.MultiPolygon(polys), poly_b; target, fix_multipoly)
        end
        #= One multipoly_a has been completely covered (and thus removed) there is no need to
        continue taking the difference =#
        isempty(polys) && break
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

