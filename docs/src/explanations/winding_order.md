# Winding order

Winding order refers specifically to the "direction" that a polygon's rings go in.  This has several different and conflicting definitions, which you can find some discussion of in the following articles:
- [GIS Stack Exchange: Order of polygon vertices?](https://gis.stackexchange.com/questions/119150/order-of-polygon-vertices-in-general-gis-clockwise-or-counterclockwise)
- [ObservableHQ winding order article](https://observablehq.com/@d3/winding-order)

GeometryOps assumes that polygon exteriors are clockwise and interiors are counterclockwise.  However, most algorithms are agnostic to winding order, and instead rely on the GeoInterface `getexterior` and `gethole` functions to distinguish holes from exteriors.  Notably, _most_ GIS implementations agree that polygons can have only one exterior but several holes.

## What other libraries do

TODO: Markdown table with a bunch of libraries/standards, their winding orders, and references.