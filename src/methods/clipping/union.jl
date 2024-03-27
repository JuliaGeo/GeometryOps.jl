# # Union Polygon Clipping
export union

"""
    union(geom_a, geom_b, [::Type{T}]; target::Type)

Return the union between two geometries as a list of geometries. Return an empty list if
none are found. The type of the list will be constrained as much as possible given the input
geometries. Furthermore, the user can provide a `taget` type as a keyword argument and a
list of target geometries found in the difference will be returned. The user can also
provide a float type 'T' that they would like the points of returned geometries to be. 
    
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
    geom_a, geom_b, ::Type{T}=Float64; target=nothing,
) where {T<:AbstractFloat}
    _union(TraitTarget(target), T, GI.trait(geom_a), geom_a, GI.trait(geom_b), geom_b)
end

_union_delay_cross_f(x) = (x, !x)
_union_delay_bounce_f(x, _) = !x
_union_step(x, _) = x ? (-1) : 1

#= This 'union' implementation returns the union of two polygons. The algorithm to determine
the union was adapted from "Efficient clipping of efficient polygons," by Greiner and
Hormann (1998). DOI: https://doi.org/10.1145/274363.274364 =#
function _union(
    ::TraitTarget{GI.PolygonTrait}, ::Type{T},
    ::GI.PolygonTrait, poly_a,
    ::GI.PolygonTrait, poly_b,
) where T
    # First, I get the exteriors of the two polygons
    ext_a = GI.getexterior(poly_a)
    ext_b = GI.getexterior(poly_b)
    # Then, I get the union of the exteriors
    a_list, b_list, a_idx_list = _build_ab_list(T, ext_a, ext_b, _union_delay_cross_f, _union_delay_bounce_f)
    polys = _trace_polynodes(T, a_list, b_list, a_idx_list, _union_step)
    n_pieces = length(polys)
    # Check if one polygon totally within other and if so, return the larger polygon.
    if n_pieces == 0 # no crossing points, determine if either poly is inside the other
        a_in_b, b_in_a = _find_non_cross_orientation(a_list, b_list, ext_a, ext_b)
        if a_in_b
            push!(polys, GI.Polygon([tuples(ext_b)]))
        elseif b_in_a
            push!(polys,  GI.Polygon([tuples(ext_a)]))
        else
            push!(polys, tuples(poly_a))
            push!(polys, tuples(poly_b))
            return polys
        end
    elseif n_pieces > 1  # extra polygons are holes (n_pieces == 1 is the desired state)
        sort!(polys, by = area, rev = true)  # sort so first element is the exterior
    end
    # the first element is the exterior, the rest are holes
    new_holes = @views (GI.getexterior(p) for p in polys[2:end])
    polys = n_pieces > 1 ? polys[1:1] : polys
    # Add holes back in for there are any
    if GI.nhole(poly_a) != 0 || GI.nhole(poly_b) != 0 || n_pieces > 1
        hole_iterator = Iterators.flatten((GI.gethole(poly_a), GI.gethole(poly_b), new_holes))
        _add_holes_to_polys!(T, polys, hole_iterator)
    end
    return polys
end


# Many type and target combos aren't implemented
function _union(
    ::TraitTarget{Target}, ::Type{T},
    trait_a::GI.AbstractTrait, geom_a,
    trait_b::GI.AbstractTrait, geom_b,
) where {Target,T}
    throw(ArgumentError("Union between $trait_a and $trait_b with target $Target isn't implemented yet."))
    return nothing
end
