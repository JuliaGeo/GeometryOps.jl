# # Covers

export covers

#=
## What is covers?

The covers function checks if a given geometry completly covers another
geometry. For this to be true, the "contained" geometry's interior and
boundaries must be covered by the "covering" geometry's interior and boundaries.
The interiors do not need to overlap.

To provide an example, consider these two lines:
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie


```

```@example cshape
covers(l1, l2)  # returns true
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
line = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])

```
"""
covers(g1, g2)::Bool = coveredby(g2, g1)
