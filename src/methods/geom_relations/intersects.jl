# # Intersection checks

export intersects

#=
## What is `intersects`?

The intersects function checks if a given geometry intersects with another
geometry, or in other words, the either the interiors or boundaries of the two
geometries intersect.

To provide an example, consider these two lines:
```@example intersects_intersection
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
f, a, p = lines(GI.getpoint(line1))
lines!(GI.getpoint(line2))
f
```
We can see that they intersect, so we expect intersects to return true, and we
can visualize the intersection point in red.
```@example intersects_intersection
GO.intersects(line1, line2)  # true
```

## Implementation

This is the GeoInterface-compatible implementation.

Given that intersects is the exact opposite of disjoint, we simply pass the two
inputs variables, swapped in order, to disjoint.
=#
"""
    intersects(geom1, geom2)::Bool

Return true if the interiors or boundaries of the two geometries interact.

`intersects` returns the exact opposite result of `disjoint`.

## Example

```jldoctest
import GeoInterface as GI, GeometryOps as GO

line1 = GI.Line([(124.584961,-12.768946), (126.738281,-17.224758)])
line2 = GI.Line([(123.354492,-15.961329), (127.22168,-14.008696)])
GO.intersects(line1, line2)

# output
true
```
"""
intersects(geom1, geom2) = !disjoint(geom1, geom2)


#= Returns true if there is at least one intersection between edges within the
two lists of edges. =#
function _line_intersects(
    edges_a::Vector{<:Edge},
    edges_b::Vector{<:Edge}
)
    # Extents.intersects(to_extent(edges_a), to_extent(edges_b)) || return false
    for edge_a in edges_a
        for edge_b in edges_b
            _line_intersects(edge_a, edge_b) && return true 
        end
    end
    return false
end

# Returns true if there is at least one intersection between two edges.
function _line_intersects(edge_a::Edge, edge_b::Edge)
    seg_val, _, _ = _intersection_point(Float64, edge_a, edge_b)
    return seg_val != line_out
end