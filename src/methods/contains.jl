# # Containment

export contains

"""
    contains(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool

Return true if the second geometry is completely contained by the first geometry.
The interiors of both geometries must intersect and, the interior and boundary of the secondary (geometry b)
must not intersect the exterior of the primary (geometry a).
`contains` returns the exact opposite result of `within`.

## Examples

```jldoctest
line = GI.LineString([[1, 1], [1, 2], [1, 3], [1, 4]])
point = Point([1, 2])
contains(line, point)
# output
true
```
"""
contains(g1, g2)::Bool = within(g2, g1)
