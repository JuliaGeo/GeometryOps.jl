# # Contains

export contains

#=
## What is contains?

The contains function checks if completly contains another geometry, or in other
words, that the second geometry is completly within the first.

To provide an example, consider these two lines:
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

l1 = GI.LineString([(0.0, 0.0), (1.0, 0.0), (0.0, 0.1)])
l2 = GI.LineString([(0.25, 0.0), (0.75, 0.0)])
f, a, p = lines(GI.getpoint(l1), color = :blue)
scatter!(GI.getpoint(l1), color = :blue)
lines!(GI.getpoint(l2), color = :orange)
scatter!(GI.getpoint(l2), color = :orange)
```
We can see that all of the points and edges of l2 are within l1, so l1 contains
l2.
```@example cshape
contains(l1, l2)  # returns true
```

## Implementation

This is the GeoInterface-compatible implementation.

Given that contains is the exact opposite of within, we simply pass the two
inputs variables, swapped in order, to within.
=#

"""
    contains(g1::AbstractGeometry, g2::AbstractGeometry)::Bool

Return true if the second geometry is completely contained by the first
geometry. The interiors of both geometries must intersect and, the interior and
boundary of the secondary (g2) must not intersect the exterior of the primary
(g1).

`contains` returns the exact opposite result of `within`.

## Examples

```jldoctest
import GeometryOps as GO, GeoInterface as GI
line = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])
point = (1, 2)

GO.contains(line, point)
# output
true
```
"""
contains(g1, g2)::Bool = within(g2, g1)
