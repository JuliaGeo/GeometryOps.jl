# #  Union Polygon Clipping
export union

"""
    union(geom_a, geom_b, ::Type{T}; target::Type{Target})

Return the union between two geometries as a list of geometries. Return an empty list if
none are found. The type of the list will be constrained as much as possible given the input
geometries. Furthermore, the user can provide a `taget` type as a keyword argument and a
list of target geometries found in the difference will be returned. The user can also
provide a float type that they would like the points of returned geometries to be. 
    
Calculates the union between two polygons.
## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
GO.union(p1, p2; target = GI.PolygonTrait)

# output
[GI.Polygon([[(6.5, 3.5), (5.0, 5.0), (0.0, 0.0), (5.0, -5.0), (6.5, -3.5), (8.0, -5.0), (13.0, 0.0), (8.0, 5.0), (6.5, 3.5)]])]
```
"""
function union(
    geom_a, geom_b, ::Type{T} = Float64; target::Type{Target} = Nothing,
) where {T <: AbstractFloat, Target <: GI.AbstractTrait}
    _union(Target, T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end

#= This 'union' implementation returns the union of two polygons. The algorithm to determine
the union was adapted from "Efficient clipping of efficient polygons," by Greiner and
Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _union(
    ::Type{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b,
) where T
    # First, I get the exteriors of the two polygons
    ext_poly_a = GI.getexterior(poly_a)
    ext_poly_b = GI.getexterior(poly_b)
    # Then, I get the union of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_poly_a, ext_poly_b)
    polys = _trace_polynodes(a_list, b_list, a_idx_list, (x, y) -> x ? (-1) : 1)
    # Check if one polygon totally within other and if so, return the larger polygon.
    if isempty(polys)
        if _point_filled_curve_orientation(a_list[1].point, ext_poly_b) == point_in
            push!(polys, GI.Polygon([ext_poly_b]))
        elseif _point_filled_curve_orientation(b_list[1].point, ext_poly_a) == point_in
            push!(polys,  GI.Polygon([ext_poly_a]))
        else
            push!(polys, poly_a)
            push!(polys, poly_b)
            return polys
        end
    else
        sort!(polys, by = area, rev = true)
        polys = [GI.Polygon([GI.getexterior(p) for p in polys])]
    end

    n_b_holes = GI.nhole(poly_b)
    if GI.nhole(poly_a) != 0 || n_b_holes != 0
        new_poly = [GI.getexterior(polys[1]); collect(GI.gethole(polys[1]))]
        current_poly = GI.Polygon([ext_poly_b])
        for (i, hole) in enumerate(Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b))))
            # Use ext_poly_b to not overcount overlapping holes in poly_a and in poly_b
            new_hole = difference(GI.Polygon([hole]), current_poly, T; target = GI.PolygonTrait)
            for h in new_hole
                push!(new_poly, GI.getexterior(h))
            end
            if i == n_b_holes
                current_poly = poly_a
            end
        end
        polys[1] = GI.Polygon(new_poly)
    end
    return polys
end

# Many type and target combos aren't implemented
function _union(
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