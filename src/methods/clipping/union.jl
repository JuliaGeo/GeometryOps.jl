# # Union Polygon Clipping
export union

"""
    union(geom_a, geom_b, [::Type{T}]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the union as a list of geometries, empty when no result exists. The list type is
constrained by the inputs; `target` selects output geometry types and `T` sets coordinate
precision.

`fix_multipoly` corrects intersecting multipolygon components before clipping. The default
union correction inherits the caller’s algorithm, manifold, and numeric type. Set
`fix_multipoly = nothing` only when the input multipolygons are valid.
## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
union_poly = GO.union(p1, p2; target = GI.PolygonTrait())
GI.coordinates.(union_poly)

# output
1-element Vector{Vector{Vector{Vector{Float64}}}}:
 [[[6.5, 3.5], [5.0, 5.0], [0.0, 0.0], [5.0, -5.0], [6.5, -3.5], [8.0, -5.0], [13.0, 0.0], [8.0, 5.0], [6.5, 3.5]]]
```
"""
function union(
    alg::FosterHormannClipping, geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...
) where {T<:AbstractFloat}
    return _union(
        alg, TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b;
        exact = True(), kwargs...,
    )
end

# fallback definitions
# if no manifold - assume planar (until we have best_manifold)
function union(
    geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...
) where {T<:AbstractFloat}
    return union(FosterHormannClipping(Planar()), geom_a, geom_b, T; target, kwargs...)
end

# if manifold but no algorithm - assume FosterHormannClipping with provided manifold.
function union(m::Manifold, geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...) where {T<:AbstractFloat}
    return union(FosterHormannClipping(m), geom_a, geom_b, T; target, kwargs...)
end

#= This 'union' implementation returns the union of two polygons. The algorithm to determine
the union was adapted from "Efficient clipping of efficient polygons," by Greiner and
Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _union(
    alg::FosterHormannClipping, ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    exact, kwargs...,
) where T
    # First, I get the exteriors of the two polygons
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then, I get the union of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(alg, T, ext_a, ext_b, _union_delay_cross_f, _union_delay_bounce_f; exact)
    polys = _trace_polynodes(alg, T, a_list, b_list, a_idx_list, _union_step, poly_a, poly_b)
    n_pieces = length(polys)
    # Check if one polygon totally within other and if so, return the larger polygon
    a_in_b, b_in_a = false, false
    if n_pieces == 0 # no crossing points, determine if either poly is inside the other
        #= Every branch here returns a piece of the input, so each is rebuilt in the
        representation `polys` is committed to rather than plain tuples. =#
        P = _fh_out_point_type(alg.manifold, poly_a, T)
        a_in_b, b_in_a = _find_non_cross_orientation(alg, a_list, b_list, ext_a, ext_b; exact)
        if a_in_b
            push!(polys, GI.Polygon([_fh_as_ring(P, ext_b, T)]))
        elseif b_in_a
            push!(polys,  GI.Polygon([_fh_as_ring(P, ext_a, T)]))
        else
            push!(polys, _fh_as_poly(P, poly_a, T))
            push!(polys, _fh_as_poly(P, poly_b, T))
            return polys
        end
    elseif n_pieces > 1
        #= extra polygons are holes (n_pieces == 1 is the desired state) and since
        holes are formed by regions exterior to both poly_a and poly_b, they can't interact
        with pre-existing holes =#
        sort!(polys, by = p -> area(alg.manifold, p, T), rev = true)  # sort by area so first element is the exterior
        # the first element is the exterior, the rest are holes
        @views append!(polys[1].geom, (GI.getexterior(p) for p in polys[2:end]))
        keepat!(polys, 1)
    end
    # Add in holes
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        _add_union_holes!(alg, polys, a_in_b, b_in_a, poly_a, poly_b; exact)
    end
    # Remove unneeded collinear points on same edge
    _remove_collinear_points!(alg, polys, [false], poly_a, poly_b)
    return polys
end

# # Helper functions for Unions with Greiner and Hormann Polygon Clipping

#= When marking the crossing status of a delayed crossing, the chain start point is crossing
when the start point is a entry point and is a bouncing point when the start point is an
exit point. The end of the chain has the opposite crossing / bouncing status. =#
_union_delay_cross_f(x) = (x, !x)

#= Delayed-bounce endpoints bounce if adjacent edges lie inside the other polygon
(`x`); otherwise they cross. =#
_union_delay_bounce_f(x, _) = !x

#= When tracing polygons, step backwards if the most recent intersection point was an entry
point, else step forwards where x is the entry/exit status. =#
_union_step(x, _) = x ? (-1) : 1

#= Add holes from two polygons to the exterior polygon formed by their union. If adding the
the holes reveals that the polygons aren't actually intersecting, return the original
polygons. =#
function _add_union_holes!(alg::FosterHormannClipping, polys, a_in_b, b_in_a, poly_a, poly_b; exact)
    P = _fh_poly_point_type(eltype(polys))
    T = _fh_float_type(P)
    if a_in_b
        _add_union_holes_contained_polys!(alg, polys, poly_a, poly_b; exact)
    elseif b_in_a
        _add_union_holes_contained_polys!(alg, polys, poly_b, poly_a; exact)
    else  # Polygons intersect, but neither is contained in the other
        n_a_holes = GI.nhole(poly_a)
        ext_poly_a = GI.Polygon(StaticArrays.SVector(GI.getexterior(poly_a)))
        ext_poly_b = GI.Polygon(StaticArrays.SVector(GI.getexterior(poly_b)))
        #= Start with poly_b when comparing with holes from poly_a and then switch to poly_a
        to compare with holes from poly_b. For current_poly, use ext_poly_b to avoid
        repeating overlapping holes in poly_a and poly_b =#
        curr_exterior_poly = n_a_holes > 0 ? ext_poly_b : ext_poly_a
        current_poly = n_a_holes > 0 ? ext_poly_b : poly_a
        # Loop over all holes in both original polygons
        for (i, ih) in enumerate(Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b))))
            ih = _fh_as_ring(P, ih, T)
            in_ext, _, _ = _line_polygon_interactions(alg.manifold, ih, curr_exterior_poly; exact, closed_line = true)
            if !in_ext
                #= if the hole isn't in the overlapping region between the two polygons, add
                the hole to the resulting polygon as we know it can't interact with any
                other holes =#
                push!(polys[1].geom, ih)
            else
                #= Subtract the other polygon from holes that intersect the overlap.
                Include its holes so shared holes within the overlap are retained. =#
                h_poly = GI.Polygon(StaticArrays.SVector(ih))
                new_holes = difference(alg, h_poly, current_poly, T; target = GI.PolygonTrait())
                append!(polys[1].geom, (GI.getexterior(new_h) for new_h in new_holes))
            end
            if i == n_a_holes
                curr_exterior_poly = ext_poly_a
                current_poly = poly_a
            end
        end
    end
    return
end

#= Add holes holes to the union of two polygons where one of the original polygons was
inside of the other. If adding the the holes reveal that the polygons aren't actually
intersecting, return the original polygons.=#
function _add_union_holes_contained_polys!(alg::FosterHormannClipping, polys, interior_poly, exterior_poly; exact)
    P = _fh_poly_point_type(eltype(polys))
    T = _fh_float_type(P)
    union_poly = polys[1]
    interior_poly = _fh_as_poly(P, interior_poly, T)
    exterior_poly = _fh_as_poly(P, exterior_poly, T)
    ext_int_ring = GI.getexterior(interior_poly)
    for (i, ih) in enumerate(GI.gethole(exterior_poly))
        poly_ih = GI.Polygon(StaticArrays.SVector(ih))
        in_ih, on_ih, out_ih = _line_polygon_interactions(alg.manifold, ext_int_ring, poly_ih; exact, closed_line = true)
        if in_ih  # at least part of interior polygon exterior is within the ith hole
            if !on_ih && !out_ih
                #= interior polygon is completely within the ith hole - polygons aren't
                touching and do not actually form a union =#
                P = _fh_poly_point_type(eltype(polys))
                polys[1] = _fh_as_poly(P, interior_poly)
                push!(polys, _fh_as_poly(P, exterior_poly))
                return polys
            else
                #= interior polygon is partially within the ith hole - area of interior
                polygon reduces the size of the hole =#
                new_holes = difference(alg, poly_ih, interior_poly, T; target = GI.PolygonTrait())
                append!(union_poly.geom, (GI.getexterior(new_h) for new_h in new_holes))
            end
        else  # none of interior polygon exterior is within the ith hole
            if !out_ih
                #= interior polygon's exterior is the same as the ith hole - polygons do
                form a union, but do not overlap so all holes stay in final polygon =#
                append!(union_poly.geom, Iterators.drop(GI.gethole(exterior_poly), i))
                append!(union_poly.geom, GI.gethole(interior_poly))
                return polys
            else
                #= interior polygon's exterior is outside of the ith hole - the interior
                polygon could either be disjoint from the hole, or contain the hole =#
                ext_int_poly = GI.Polygon(StaticArrays.SVector(ext_int_ring))
                in_int, _, _ = _line_polygon_interactions(alg.manifold, ih, ext_int_poly; exact, closed_line = true)
                if in_int
                    #= interior polygon contains the hole - overlapping holes between the
                    interior and exterior polygons will be added =#
                    for jh in GI.gethole(interior_poly)
                        poly_jh = GI.Polygon(StaticArrays.SVector(jh))
                        if intersects(alg.manifold, poly_ih, poly_jh)
                            new_holes = intersection(alg, poly_ih, poly_jh, T; target = GI.PolygonTrait())
                            append!(union_poly.geom, (GI.getexterior(new_h) for new_h in new_holes))
                        end
                    end
                else
                    #= interior polygon and the exterior polygon are disjoint - add the ith
                    hole as it is not covered by the interior polygon =#
                    push!(union_poly.geom, ih)
                end
            end
        end
    end
    return
end

#= Include all components of `multipoly_b`, merging those that intersect `poly_a`.
Correct the multipolygon unless `fix_multipoly = nothing`. =#
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent repeated regions in the output
        multipoly_b = fix_multipoly(multipoly_b)
    end
    P = _fh_out_point_type(alg.manifold, poly_a, T)
    polys = [_fh_as_poly(P, poly_a, T)]
    for poly_b in GI.getpolygon(multipoly_b)
        if intersects(alg.manifold, polys[1], poly_b)
            # If polygons intersect and form a new polygon, swap out polygon
            new_polys = union(alg, polys[1], poly_b, T; target)
            if length(new_polys) > 1 # case where they intersect by just one point
                push!(polys, _fh_as_poly(P, poly_b, T))  # add poly_b to list
            else
                polys[1] = new_polys[1]
            end
        else
            # If they don't intersect, poly_b is now a part of the union as its own polygon
            push!(polys, _fh_as_poly(P, poly_b, T))
        end
    end
    return polys
end

#= Preserve operand order so mixed representations do not round-trip passthrough vertices.
The multipolygon implementation accumulates its original components in their own chart. =#
_union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T = union(alg, multipoly_a, GI.MultiPolygon([poly_b]), T; target, kwargs...)

#= Merge intersecting components of both multipolygons. Apply the supplied correction
unless `fix_multipoly = nothing`. =#
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent repeated regions in the output
        multipoly_a = fix_multipoly(multipoly_a)
        multipoly_b = fix_multipoly(multipoly_b)
        fix_multipoly = nothing
    end
    multipolys = multipoly_b
    polys = [_fh_as_poly(_fh_out_point_type(alg.manifold, multipoly_a, T), p, T) for p in GI.getpolygon(multipoly_b)]
    for poly_a in GI.getpolygon(multipoly_a)
        polys = union(alg, poly_a, multipolys, T; target, fix_multipoly)
        multipolys = _fh_multipolygon(polys)
    end
    return polys
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T
    return _fh_multipolygon(union(alg, poly_a, poly_b, T; target = GI.PolygonTrait(), kwargs...))
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T),
    kwargs...,
) where T 
    res = union(alg, poly_a, multipoly_b, T; target = GI.PolygonTrait(), fix_multipoly, kwargs...)
    if !isnothing(fix_multipoly)
        return fix_multipoly(_fh_multipolygon(res))
    else
        return _fh_multipolygon(res)
    end
end
# this is the opposite of the above
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T),
    kwargs...,
) where T 
    return _fh_multipolygon(union(alg, multipoly_a, poly_b, T;
        target = GI.PolygonTrait(), fix_multipoly, kwargs...))
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    trait_a::GI.MultiPolygonTrait, multipoly_a,
    trait_b::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(alg, T), kwargs...,
) where T
    return _fh_multipolygon(_union(alg, TraitTarget{GI.PolygonTrait}(), T, trait_a, multipoly_a, trait_b, multipoly_b; fix_multipoly, kwargs...))
end
# Many type and target combos aren't implemented
function _union(
    alg::GeometryOpsCore.Algorithm, target::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b;
    kwargs...
) where {Target,T}
    throw(ArgumentError("Union between $trait_a and $trait_b with target $Target isn't implemented yet."))
    return nothing
end
