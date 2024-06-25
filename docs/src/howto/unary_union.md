# Unary union

```@example unary
using Makie, GeoMakie
import GeometryOps as GO, GeoInterface as GI
```

Furst, we get the data.  These are watersheds in Vancouver Island, Canada - the very same used in the GEOS benchmarks.

Since the file is zipped, we have to unzip it, which is why this code is a bit longer than it otherwise might be.

```@example unary
wkt_gz_file = download("https://rawcdn.githack.com/pramsey/geos-performance/b54d92a678e2174059d1b0ff233e275e4bd02084/data/watersheds.wkt.gz", joinpath(tempdir(), "watersheds.wkt.gz"))
using GZip
handle = GZip.open(wkt_gz_file)

using WellKnownGeometry, GeoFormatTypes
wkt = GeoFormatTypes.WellKnownText.((GeoFormatTypes.Geom(),), eachline(handle))
close(handle)

geoms = GO.tuples(wkt)

plot(geoms; color = 1:length(geoms), axis = (; aspect = DataAspect()))
```

Now that we have the geometries, we reduce over the vector, performing unions along the way.

```@example unary
@time final_multipoly = reduce(
    (x, y) -> GO.union(x, y; target = GI.MultiPolygonTrait, fix = GO.UnionIntersectingPolygons()), 
    GO.fix(geoms)
)
```


```@example unary
fixed_geoms = GO.buffer(geoms, 0)
@time final_multipoly = reduce(
    (x, y) -> GO.union(x, y; target = GI.MultiPolygonTrait, fix = GO.UnionIntersectingPolygons()), 
    fixed_geoms
)
```
```@example unary
@time GO.fix(GI.MultiPolygon(fixed_geoms); corrections = (GO.UnionIntersectingPolygons(),))
```
