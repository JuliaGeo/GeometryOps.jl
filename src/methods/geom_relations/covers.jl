# # Covers

export covers

#=
## What is covers?

The covers function checks if a given geometry completely covers another
geometry. For this to be true, the "contained" geometry's interior and
boundaries must be covered by the "covering" geometry's interior and boundaries.
The interiors do not need to overlap.

To provide an example, consider these two lines:
```@example covers
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

p1 = (0.0, 0.0)
p2 = (1.0, 1.0)
l1 = GI.Line([p1, p2])

f, a, p = lines(GI.getpoint(l1))
scatter!(p1, color = :red)
f
```

```@example covers
GO.covers(l1, p1)  # returns true
GO.covers(p1, l1)  # returns false
```

## Implementation

This is the GeoInterface-compatible implementation.

Given that covers is the exact opposite of coveredby, we simply pass the two
inputs variables, swapped in order, to coveredby.
=#

"""
    covers(g1::AbstractGeometry, g2::AbstractGeometry)::Bool

Return true if the first geometry is completely covers the second geometry, 
The exterior and boundary of the second geometry must not be outside of the
interior and boundary of the first geometry. However, the interiors need not
intersect.

`covers` returns the exact opposite result of `coveredby`.

## Examples

```jldoctest
import GeometryOps as GO, GeoInterface as GI
l1 = GI.LineString([(1.0, 1.0), (1.0, 2.0), (1.0, 3.0), (1.0, 4.0)])
l2 = GI.LineString([(1.0, 1.0), (1.0, 2.0)])

GO.covers(l1, l2)
# output
true
```
"""
covers(g1, g2)::Bool = GeometryOps.coveredby(g2, g1)

"""
    covers(g1)

Return a function that checks if its input covers `g1`.
This is equivalent to `x -> covers(x, g1)`.
"""
covers(g1) = Base.Fix2(covers, g1)