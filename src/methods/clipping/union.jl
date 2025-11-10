# # Union Polygon Clipping
export union

"""
    union(geom_a, geom_b, [::Type{T}]; target::Type, fix_multipoly = UnionIntersectingPolygons())

Return the union between two geometries as a list of geometries. Return an empty list if
none are found. The type of the list will be constrained as much as possible given the input
geometries. Furthermore, the user can provide a `taget` type as a keyword argument and a
list of target geometries found in the difference will be returned. The user can also
provide a float type 'T' that they would like the points of returned geometries to be. If
the user is taking a intersection involving one or more multipolygons, and the multipolygon
might be comprised of polygons that intersect, if `fix_multipoly` is set to an
`IntersectingPolygons` correction (the default is `UnionIntersectingPolygons()`), then the
needed multipolygons will be fixed to be valid before performing the intersection to ensure
a correct answer. Only set `fix_multipoly` to false if you know that the multipolygons are
valid, as it will avoid unneeded computation. 
    
Calculates the union between two polygons.
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
    return union(FosterHormannClipping(Planar()), geom_a, geom_b; target, kwargs...)
end

# if manifold but no algorithm - assume FosterHormannClipping with provided manifold.
function union(m::Manifold, geom_a, geom_b, ::Type{T}=Float64; target=nothing, kwargs...) where {T<:AbstractFloat}
    return union(FosterHormannClipping(m), geom_a, geom_b; target, kwargs...)
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

    crs = mutual_crs(poly_a, poly_b)
    # Then, I get the union of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(alg, T, ext_a, ext_b, _union_delay_cross_f, _union_delay_bounce_f; exact)
    polys = _trace_polynodes(alg, T, a_list, b_list, a_idx_list, _union_step, poly_a, poly_b; crs)
    n_pieces = length(polys)
    # Check if one polygon totally within other and if so, return the larger polygon
    a_in_b, b_in_a = false, false
    if n_pieces == 0 # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(alg, a_list, b_list, ext_a, ext_b; exact)
        if a_in_b
            push!(polys, GI.Polygon([_linearring(tuples(ext_b; crs))]; crs))
        elseif b_in_a
            push!(polys,  GI.Polygon([_linearring(tuples(ext_a; crs))]; crs))
        else
            push!(polys, tuples(poly_a; crs))
            push!(polys, tuples(poly_b; crs))
            return polys
        end
    elseif n_pieces > 1
        #= extra polygons are holes (n_pieces == 1 is the desired state) and since
        holes are formed by regions exterior to both poly_a and poly_b, they can't interact
        with pre-existing holes =#
        sort!(polys, by = area, rev = true)  # sort by area so first element is the exterior
        # the first element is the exterior, the rest are holes
        @views append!(polys[1].geom, (GI.getexterior(p) for p in polys[2:end]))
        keepat!(polys, 1)
    end
    # Add in holes
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        _add_union_holes!(alg, polys, a_in_b, b_in_a, poly_a, poly_b; exact, crs)
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

#= When marking the crossing status of a delayed bouncing, the chain start and end points
are bouncing if the current polygon's adjacent edges are within the non-tracing polygon. If
the edges are outside then the chain endpoints are marked as crossing. x is a boolean
representing if the edges are inside or outside of the polygon. =#
_union_delay_bounce_f(x, _) = !x

#= When tracing polygons, step backwards if the most recent intersection point was an entry
point, else step forwards where x is the entry/exit status. =#
_union_step(x, _) = x ? (-1) : 1

#= Add holes from two polygons to the exterior polygon formed by their union. If adding the
the holes reveals that the polygons aren't actually intersecting, return the original
polygons. =#
function _add_union_holes!(alg::FosterHormannClipping, polys, a_in_b, b_in_a, poly_a, poly_b; exact, crs)
    if a_in_b
        _add_union_holes_contained_polys!(alg, polys, poly_a, poly_b; exact, crs)
    elseif b_in_a
        _add_union_holes_contained_polys!(alg, polys, poly_b, poly_a; exact, crs)
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
            ih = _linearring(ih)
            in_ext, _, _ = _line_polygon_interactions(#=TODO: alg.manifold=# ih, curr_exterior_poly; exact, closed_line = true)
            if !in_ext
                #= if the hole isn't in the overlapping region between the two polygons, add
                the hole to the resulting polygon as we know it can't interact with any
                other holes =#
                push!(polys[1].geom, ih)
            else
                #= if the hole is at least partially in the overlapping region, take the
                difference of the hole from the polygon it didn't originate from - note that
                when current_poly is poly_a this includes poly_a holes so overlapping holes
                between poly_a and poly_b within the overlap are added, in addition to all
                holes in non-overlapping regions =#
                h_poly = GI.Polygon(StaticArrays.SVector(ih))
                new_holes = difference(alg, h_poly, current_poly; target = GI.PolygonTrait())
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
function _add_union_holes_contained_polys!(alg::FosterHormannClipping, polys, interior_poly, exterior_poly; exact, crs)
    union_poly = polys[1]
    ext_int_ring = GI.getexterior(interior_poly)
    for (i, ih) in enumerate(GI.gethole(exterior_poly))
        poly_ih = GI.Polygon(StaticArrays.SVector(ih))
        in_ih, on_ih, out_ih = _line_polygon_interactions(#=TODO: alg.manifold=# ext_int_ring, poly_ih; exact, closed_line = true)
        if in_ih  # at least part of interior polygon exterior is within the ith hole
            if !on_ih && !out_ih
                #= interior polygon is completely within the ith hole - polygons aren't
                touching and do not actually form a union =#
                polys[1] = tuples(interior_poly; crs)
                push!(polys, tuples(exterior_poly; crs))
                return polys
            else
                #= interior polygon is partially within the ith hole - area of interior
                polygon reduces the size of the hole =#
                new_holes = difference(alg, poly_ih, interior_poly; target = GI.PolygonTrait())
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
                in_int, _, _ = _line_polygon_interactions(#=TODO: alg.manifold=# ih, ext_int_poly; exact, closed_line = true)
                if in_int
                    #= interior polygon contains the hole - overlapping holes between the
                    interior and exterior polygons will be added =#
                    for jh in GI.gethole(interior_poly)
                        poly_jh = GI.Polygon(StaticArrays.SVector(jh))
                        if intersects(poly_ih, poly_jh)
                            new_holes = intersection(#=TODO: alg, =#poly_ih, poly_jh; target = GI.PolygonTrait())
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

#= Polygon with multipolygon union - note that all sub-polygons of `multipoly_b` will be
included, unioning these sub-polygons with `poly_a` where they intersect. Unless specified
with `fix_multipoly = nothing`, `multipolygon_b` will be validated using the given (default
is `UnionIntersectingPolygons()`) correction. =#
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent repeated regions in the output
        multipoly_b = fix_multipoly(multipoly_b)
    end
    crs = mutual_crs(poly_a, multipoly_b)
    polys = [tuples(poly_a, T; crs)]
    for poly_b in GI.getpolygon(multipoly_b)
        if intersects(#=TODO: alg.manifold, =#polys[1], poly_b)
            # If polygons intersect and form a new polygon, swap out polygon
            new_polys = union(alg, polys[1], poly_b; target)
            if length(new_polys) > 1 # case where they intersect by just one point
                push!(polys, tuples(poly_b, T; crs))  # add poly_b to list
            else
                polys[1] = new_polys[1]
            end
        else
            # If they don't intersect, poly_b is now a part of the union as its own polygon
            push!(polys, tuples(poly_b, T; crs))
        end
    end
    return polys
end

#= Multipolygon with polygon union is equivalent to taking the union of the polygon with the
multipolygon and thus simply switches the order of operations and calls the above method. =#
_union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T = union(alg, poly_b, multipoly_a; target, kwargs...)

#= Multipolygon with multipolygon union - note that all of the sub-polygons of `multipoly_a`
and the sub-polygons of `multipoly_b` are included and combined together where there are
intersections. Unless specified with `fix_multipoly = nothing`, `multipolygon_b` will be
validated using the given (default is `UnionIntersectingPolygons()`) correction. =#
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    if !isnothing(fix_multipoly) # Fix multipoly_b to prevent repeated regions in the output
        multipoly_b = fix_multipoly(multipoly_b)
        fix_multipoly = nothing
    end
    multipolys = multipoly_b
    local polys
    for poly_a in GI.getpolygon(multipoly_a)
        polys = union(alg, poly_a, multipolys; target, fix_multipoly)
        multipolys = GI.MultiPolygon(polys)
    end
    return polys
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b;
    kwargs...,
) where T
    return UnionIntersectingPolygons(#=TODO: alg.manifold=#)(GI.MultiPolygon([poly_a, poly_b]))
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(),
    kwargs...,
) where T 
    res = union(alg, multipoly_b, poly_a; target = GI.PolygonTrait(), kwargs...)
    if !isnothing(fix_multipoly)
        return fix_multipoly(GI.MultiPolygon(res))
    else
        return GI.MultiPolygon(res)
    end
end
# this is the opposite of the above
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    ::GI.MultiPolygonTrait, multipoly_a,
    ::GI.PolygonTrait, poly_b;
    fix_multipoly = UnionIntersectingPolygons(),
    kwargs...,
) where T 
    union(alg, poly_b, multipoly_a; target, fix_multipoly, kwargs...)
end
function _union(
    alg::FosterHormannClipping, target::TraitTarget{GI.MultiPolygonTrait}, ::Type{T},
    trait_a::GI.MultiPolygonTrait, multipoly_a,
    trait_b::GI.MultiPolygonTrait, multipoly_b;
    fix_multipoly = UnionIntersectingPolygons(), kwargs...,
) where T
    return GI.MultiPolygon(_union(alg, TraitTarget{GI.PolygonTrait}(), T, trait_a, multipoly_a, trait_b, multipoly_b; fix_multipoly, kwargs...))
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
