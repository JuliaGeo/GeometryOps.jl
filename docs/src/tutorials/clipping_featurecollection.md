# Clipping a FeatureCollection by a polygon

Often, one might want to clip a feature collection (maybe from Shapefile or GeoJSON)

```@example clipping_featurecollection
function clip_or_empty(polygon, clipper)
    #=
    result = GO.intersection(polygon, clipper; target = GI.PolygonTrait())
    if isempty(result)
        null_point = GI.is3d(polygon) ? (GI.ismeasured(polygon) ? (NaN, NaN, NaN, NaN) : (NaN, NaN, NaN)) : (NaN, NaN)
        contents = GI.LinearRing.([[null_point, null_point, null_point]])
        return GI.Polygon{GI.is3d(polygon),GI.ismeasured(polygon),typeof(contents),Nothing, typeof(GI.crs(polygon))}(contents, nothing, GI.crs(polygon))
    else
        return GI.MultiPolygon(result; crs = GI.crs(polygon))
    end
    =#
    return GO.intersection(GO.GEOS(), polygon, clipper)
end
```
First, let's load our data:
```@example clipping_featurecollection
df = nothing # DataFrame(Shapefile.Table(...))
```
and plot it:
```@example clipping_featurecollection
f, a, p = poly(df.geometry)
```
Now, we can define some polygon in that space, that we want to use to clip all geometries by!
```@example clipping_featurecollection
clipping_poly = GI.Polygon([[(880_000, 990_000), (910_000, 990_000), (910_000, 1030_000), (880_000, 1030_000), (880_000, 990_000)]])
poly!(a, clipping_poly; color = Makie.Cycled(2))
f
```
Finally, we clip, and show the output:
```@example clipping_featurecollection
clipped_geoms = clip_or_empty.(df.geometry, (clipping_poly,))
```
```@example clipping_featurecollection
poly!(a, clipped_geoms; color = Makie.Cycled(3), strokewidth = 0.75, strokecolor = :forestgreen)
f
```
