# Geodesic paths

Geodesic paths are paths computed on an ellipsoid, as opposed to a plane.  

```@example geodesic
import GeometryOps as GO, GeoInterface as GI
using CairoMakie, GeoMakie


IAH = (-95.358421, 29.749907)
AMS = (4.897070, 52.377956)


fig, ga, _cp = lines(GeoMakie.coastlines(); axis = (; type = GeoAxis))
lines!(ga, GO.segmentize(GO.GeodesicSegments(; max_distance = 100_000), GI.LineString([IAH, AMS])); color = Makie.wong_colors()[2])
fig
```