# # Difference Polygon Clipping
export difference


"""
    difference(geom_a, geom_b, [T::Type]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the difference as a list of geometries, empty when no result exists. The list type is
constrained by the inputs; `target` selects output geometry types and `T` sets coordinate
precision.

`fix_multipoly` corrects intersecting multipolygon components before clipping. The default
union correction inherits the caller’s algorithm, manifold, and numeric type. Set
`fix_multipoly = nothing` only when the input multipolygons are valid.

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
    alg::FosterHormannClipping, geom_a, geom_b, ::Type{T} = Float64; target=nothing, kwargs...,
) where {T<:AbstractFloat}
    return _difference(
        alg, TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b;
        exact = True(), kwargs...,
    )
end
# fallback definitions
difference(geom_a, geom_b, ::Type{T} = Float64; target=nothing, kwargs...) where T = difference(FosterHormannClipping(Planar()), geom_a, geom_b, T; target, kwargs...)
# if manifold but no algorithm - assume FosterHormannClipping with provided manifold.
difference(m::Manifold, geom_a, geom_b, ::Type{T} = Float64; target=nothing, kwargs...) where T = difference(FosterHormannClipping(m), geom_a, geom_b, T; target, kwargs...)

#= The 'difference' function returns the difference of two polygons as a list of polygons.
The algorithm to determine the difference was adapted from "Efficient clipping of efficient
polygons," by Greiner and Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _difference(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    exact, kwargs...
) where T
    # Get the exterior of the polygons
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Find the difference of the exterior of the polygons
    a_list, b_list, a_idx_list = _build_ab_list(alg, T, ext_a, ext_b, _diff_delay_cross_f, _diff_delay_bounce_f; exact)
    polys = _trace_polynodes(alg, T, a_list, b_list, a_idx_list, _diff_step, poly_a, poly_b)
    # if no crossing points, determine if either poly is inside of the other
    if isempty(polys)
        #= Both answers here are rebuilt from the inputs, so they take the representation
        `polys` is committed to rather than plain tuples. =#
        P = _fh_out_point_type(alg.manifold, poly_a, T)
        a_in_b, b_in_a = _find_non_cross_orientation(alg.manifold, a_list, b_list, ext_a, ext_b; exact)
        # add case for if they polygons are the same (all intersection points!)
        # add a find_first check to find first non-inter poly!
        if b_in_a && !a_in_b  # b in a and can't be the same polygon
            poly_a_b_hole = GI.Polygon([_fh_as_ring(P, ext_a, T), _fh_as_ring(P, ext_b, T)])
            push!(polys, poly_a_b_hole)
        elseif !b_in_a && !a_in_b # polygons don't intersect
            push!(polys, _fh_as_poly(P, poly_a, T))
            return polys
        end
    end
    remove_idx = falses(length(polys))
    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0
        _add_holes_to_polys!(alg, T, polys, GI.gethole(poly_a), remove_idx; exact)
    end
    if GI.nhole(poly_b) != 0
        for hole in GI.gethole(poly_b)
            hole_poly = GI.Polygon(StaticArrays.SVector(hole))
            new_polys = intersection(alg, hole_poly, poly_a, T; target = GI.PolygonTrait)
            if length(new_polys) > 0
                append!(polys, new_polys)
            end
        end
    end
    # Remove unneeded collinear points on same edge
    _remove_collinear_points!(alg, polys, remove_idx, poly_a, poly_b)
    return polys
end

# # Helper functions for Differences with Greiner and Hormann Polygon Clipping

#= When marking the crossing status of a delayed crossing, the chain start point is crossing
when the start point is a entry point and is a bouncing point when the start point is an
exit point. The end of the chain has the opposite crossing / bouncing status. =#
_diff_delay_cross_f(x) = (x, !x)
#= Delayed-bounce endpoints cross when `x ⊻ y`: `x` means inside the other polygon,
and `y` means tracing `a_list`. Otherwise they bounce. =#
_diff_delay_bounce_f(x, y) = x ⊻ y
#= Step forward for entry on `b_list` or exit on `a_list`; otherwise step backward.
`x` is entry status and `y` is true for `a_list`. =#
_diff_step(x, y) = (x ⊻ y) ? 1 : (-1)

# Subtract every component of `multipoly_b` from `poly_a`.
function _difference(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    # Although redundant subtraction is set-theoretically harmless, revisiting a
    # spherical boundary after rounding its crossings can create duplicate trace nodes.
    isnothing(fix_multipoly) || (multipoly_b = fix_multipoly(multipoly_b))
    polys = [_fh_as_poly(_fh_out_point_type(alg.manifold, poly_a, T), poly_a, T)]
    for poly_b in GI.getpolygon(multipoly_b)
        isempty(polys) && break
        polys = mapreduce(p -> difference(alg, p, poly_b, T; target), append!, polys)
    end
    return polys
end

#= Subtract `poly_b` from every component of `multipoly_a`. Apply `fix_multipoly`
unless it is `nothing`. =#
function _difference(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_a to prevent returning an invalid multipolygon
        multipoly_a = fix_multipoly(multipoly_a)
    end
    polys = Vector{_get_poly_type(T, _fh_out_point_type(alg.manifold, multipoly_a, T))}()
    sizehint!(polys, GI.npolygon(multipoly_a))
    for poly_a in GI.getpolygon(multipoly_a)
        append!(polys, difference(alg, poly_a, poly_b, T; target))
    end
    return polys
end

#= Subtract all components of `multipoly_b` from `multipoly_a`. Correct both
multipolygons unless `fix_multipoly = nothing`. =#
function _difference(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_a to prevent returning an invalid multipolygon
        multipoly_a = fix_multipoly(multipoly_a)
        multipoly_b = fix_multipoly(multipoly_b)
        fix_multipoly = nothing
    end
    polys = [_fh_as_poly(_fh_out_point_type(alg.manifold, multipoly_a, T), p, T) for p in GI.getpolygon(multipoly_a)]
    for (i, poly_b) in enumerate(GI.getpolygon(multipoly_b))
        # Subtract each component from the remaining pieces in `polys`.
        polys = if i == 1
            difference(alg, multipoly_a, poly_b, T; target, fix_multipoly)
        else
            difference(alg, _fh_multipolygon(polys), poly_b, T; target, fix_multipoly)
        end
        #= One multipoly_a has been completely covered (and thus removed) there is no need to
        continue taking the difference =#
        isempty(polys) && break
    end
    return polys
end
function _difference(
    alg::FosterHormannClipping, ::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    trait_a::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, polylike_a,
    trait_b::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, polylike_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...
) where T
    polys = _difference(alg, TraitTarget(GI.PolygonTrait()), T, trait_a, polylike_a, trait_b, polylike_b; fix_multipoly, kwargs...)
    if isnothing(fix_multipoly)
        return _fh_multipolygon(polys)
    else
        return fix_multipoly(_fh_multipolygon(polys))
    end
end
# Many type and target combos aren't implemented
function _difference(
    alg::GeometryOpsCore.Algorithm, target::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
    kw...
) where {Target, T}
    @assert(
        false,
        "Difference between $trait_a and $trait_b with target $Target and algorithm $alg isn't implemented yet.",
    )
    return nothing
end

