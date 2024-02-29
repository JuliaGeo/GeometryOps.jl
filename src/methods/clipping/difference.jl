# #  Difference Polygon Clipping
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
diff_poly = GO.difference(poly1, poly2; target = GI.PolygonTrait)
GI.coordinates.(diff_poly)

# output
1-element Vector{Vector{Vector{Vector{Float64}}}}:
 [[[6.5, 3.5], [5.0, 5.0], [0.0, 0.0], [5.0, -5.0], [6.5, -3.5], [3.0, 0.0], [6.5, 3.5]]]
```
"""
function difference(
    geom_a, geom_b, ::Type{T} = Float64; target::Type{Target} = Nothing,
) where {T <: AbstractFloat, Target <: Union{Nothing, GI.AbstractTrait}}
    return _difference(Target, T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end

#= The 'difference' function returns the difference of two polygons as a list of polygons.
The algorithm to determine the difference was adapted from "Efficient clipping of efficient
polygons," by Greiner and Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _difference(
    ::Type{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b,
) where T
    # Get the exterior of the polygons
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Find the difference of the exterior of the polygons
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, (x, y) -> (x ⊻ y) ? 1 : (-1))
    # if no crossing points, determine if either poly is inside of the other
    if isempty(polys)
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b)
        # add case for if they polygons are the same (all intersection points!)
        # add a find_first check to find first non-inter poly!
        if b_in_a && !a_in_b  # b in a and can't be the same polygon
            share_edge_warn(a_list, "Edge case: polygons share edge but one is hole of the other.")  # will get taken care of with "glued edges"
            poly_a_b_hole = GI.Polygon([tuples(ext_a), tuples(ext_b)])
            push!(polys, poly_a_b_hole)
        elseif !b_in_a && !a_in_b # polygons don't intersect
            push!(polys, tuples(poly_a))
            return polys
        end
    end

    # If the original polygons had holes, take that into account.
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0
        _add_holes_to_polys!(T, polys, GI.gethole(poly_a))
        for hole in GI.gethole(poly_b)
            new_polys = intersection(GI.Polygon([hole]), poly_a, T; target = GI.PolygonTrait)
            if length(new_polys) > 0
                append!(polys, new_polys)
            end
        end
    end
    return polys
end

# Many type and target combos aren't implemented
function _difference(
    ::Type{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
) where {Target, T}
    @assert(
        false,
        "Difference between $trait_a and $trait_b with target $Target isn't implemented yet.",
    )
    return nothing
end

