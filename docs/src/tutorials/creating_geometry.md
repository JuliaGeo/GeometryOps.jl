# Creating Geometry

````@example creating_geometry
import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import CoordinateTransformations
import Proj
using CairoMakie
using GeoMakie
using GeoJSON
````

The first thing we need to do is decide which Coordinate Reference System (CRS) we will be working in. Here, we start with the most common geographic CRS (i.e. coordiantes of latitude and longitude), [WGS84](https://epsg.io/4326).

````@example creating_geometry
crs = GFT.EPSG(4326)
````

Let's start by making a single `Point`.

````@example creating_geometry
point = GI.Point(0, 0)
````

Now, let's plot our point.

````@example creating_geometry
fig, ax, plt = plot(point)
````

Let's create a set of points, and have a bit more fun with plotting.

````@example creating_geometry
x = [-5, 0, 5, 0]
y = [0, -5, 0, 5]
points = GI.Point.(zip(x,y))
plot!(ax, points; marker = '✈', markersize = 30)
fig
````

Points can be combined into a single `MultiPoint` geometry. This time let's include information on CRS with the geometry making it a geospatial geometry. All that's needed is to include `; crs = crs` as a keyword argument when constucting the geometry.

!!! note
    It is good practice to only include CRS information with the highest-level geometry. Not doing so can bloat the memory footprint of the geometry. CRS information _can_ be included at the individual `Point` level but is discouraged.

````@example creating_geometry
x = [-5, -5, 5, 5]
y = [-5, 5, 5, -5]
multipoint = GI.MultiPoint(GI.Point.(zip(x, y)); crs)
plot!(ax, multipoint; marker = '☁', markersize = 30)
fig
````

Let's create a line between two points.

````@example creating_geometry
p1 = GI.Point.(-5, 0)
p2 = GI.Point.(5, 0)
line = GI.LineString([p1,p2]; crs)
plot!(ax, line)
fig
````

Now, let's create a line connecting multiple points (i.e. a `LineString`).
This time we get a bit more fancy with point creation.

````@example creating_geometry
r = 2;
k = 10;
ϴ = 0:0.01:2pi
x = r .* (k + 1) .* cos.(ϴ) .- r .* cos.((k + 1) .* ϴ)
y = r .* (k + 1) .* sin.(ϴ) .- r .* sin.((k + 1) .* ϴ)
lines = GI.LineString(GI.Point.(zip(x,y)); crs)
plot!(ax, lines; linewidth = 3)
fig
````

We can also create a single `LinearRing` trait, the building block of a polygon.
A `LinearRing` is simply a `LineString` with the same beginning and endpoint,
i.e., an arbitrary closed shape composed of point pairs.

A `LinearRing` is composed of a series of points listed in clockwise order (i.e., winding order).
I always think of a polygon as filled to the right of the lines as one progresses
from point `n` to point `n+1`.

````@example creating_geometry
ring1 = GI.LinearRing(GI.getpoint(lines))
````

Now, let's make the `LinearRing` into a `Polygon`.

````@example creating_geometry
polygon1 = GI.Polygon([ring1]; crs)
````

Now, we can use GeometryOperations and CoordinateTransformations to shift `polygon1`
vertically up, to avoid plotting over our earlier results.


:::tabs

== CoordinateTransformations

````@example creating_geometry
xoffset = 0.
yoffset = 50.
f = CoordinateTransformations.Translation(xoffset, yoffset)
polygon1 = GO.transform(f, polygon1)
plot!(polygon1)
fig
````

== Direct manipulation

```julia
polygon1 = GO.transform(polygon1) do point
    # Here, you must return an object that has the 
    # `GI.PointTrait` trait.  GeoInterface interprets
    # tuples of Float64 as Points by default, so we return
    # that.
    return (GI.x(point), GI.y(point) + 50)
end
plot!(polygon1)
fig
```
```@example creating_geometry
fig # hide
```
:::

Polygons can contain "holes". The first `LinearRing` in a polygon is the exterior, and all 
subsequent `LinearRing`s are treated as holes in the leading `LinearRing`.

`GeoInterface` offers the `GI.getexterior(poly)` and `GI.gethole(poly)` methods to get the 
exterior ring and an iterable of holes, respectively.

!!! note
    Some packages always consider the secondary `LinearRings` holes, others look at the winding
    order, where the polygons are filled inward if they have a [clockwise](@ref GeometryOps.isclockwise) winding order and 
    outward if they have a counterclockwise winding order.
    
    Hopefully, these are details that you'll never have to deal with.  But it is good to know.

````@example creating_geometry
hole = GI.LinearRing(GI.getpoint(multipoint))
polygon1 = GI.Polygon([ring1, hole]; crs)
````

Shift `polygon1` to the right, to avoid plotting over our earlier results.

:::tabs

== CoordinateTransformations

````@example creating_geometry
xoffset = 50.
yoffset = 0.
f = CoordinateTransformations.Translation(xoffset, yoffset)
polygon1 = GO.transform(f, polygon1)
plot!(polygon1)
fig
````

== Direct manipulation

```julia
polygon1 = GO.transform(polygon1) do point
    return (GI.x(point) + 50, GI.y(point))
end
plot!(polygon1)
fig
```
```@example creating_geometry
fig # hide
```
:::

`Polygon`s can also be grouped together as a `MultiPolygon`.

````@example creating_geometry
r = 5
x = cos.(reverse(ϴ)) .* r .+ xoffset
y = sin.(reverse(ϴ)) .* r .+ yoffset
ring2 =  GI.LinearRing(GI.Point.(zip(x,y)))
polygon2 = GI.Polygon([ring2])
multipolygon = GI.MultiPolygon([polygon1, polygon2]; crs)
````

Shift `multipolygon` up, to avoid plotting over our earlier results.

````@example creating_geometry
xoffset = 0.
yoffset = 50.
f = CoordinateTransformations.Translation(xoffset, yoffset)
multipolygon = GO.transform(f, multipolygon)
plot!(multipolygon)
fig
````

Great, now we can make `Points`, `MultiPoints`, `Lines`, `LineStrings`, `Polygons` (with holes), and `MultiPolygons` and modify them using [`CoordinateTransformations`] and [`GeometryOps`].

But where does the `crs` information come in? To show this, we need to use `GeoMakie` that can interpret the `crs` information that we've included with our geometries.

Now specify the source and destination projections for our map. Remember that the very first thing we did was select our source coordinate system.

````@example creating_geometry
source = crs;
dest = "+proj=natearth2" #see [https://proj.org/en/9.4/operations/projections/natearth2.html]
````

Open the Natural Earth continental outlines, which are available from https://www.naturalearthdata.com/, and are bundled with GeoMakie.

````@example creating_geometry
land_path = GeoMakie.assetpath("ne_110m_land.geojson")
````

Read the land polygons into a `GeoJSON.FeatureCollection`.

````@example creating_geometry
land_geo = GeoJSON.read(read(land_path, String))
````

create a figure with a `GeoAxis` from GeoMakie, that can handle the projections between CRS.

````@example creating_geometry
fig = Figure(size=(1000, 500));
ga = GeoAxis(
    fig[1, 1];
    source=crs, # `source` and `dest` set the CRS
    dest=dest,
    xticklabelsvisible = false,
    yticklabelsvisible = false,
);
nothing #hide
````

Plot `land` for context.

````@example creating_geometry
poly!(ga, land_geo, color=:black)
fig
````

Now let's plot a `Polygon` like before, but now on coordinate reference system (CRS) that is different from our data

````@example creating_geometry
plot!(multipolygon; color = :green)
fig
````

OK, now what if we want to plot geomitires that are in a diffent Coordinate Reference System (CRS?)

Select out new source `CRS`: WGS 84 / UTM zone 10N coodinate system [EPSG:32610](https://epsg.io/32610)
````@example creating_geometry
crs2 = GFT.EPSG(32610)
````

Create a polygon in the new coordinate systen (we're working in meters no, not lat/lon)
````@example creating_geometry
r = 1000000;
ϴ = 0:0.01:2pi
x = r .* cos.(ϴ).^3 .+ 500000
y = r .* sin.(ϴ) .^ 3 .+5000000
````

Now create a `LinearRing` from `Points`

````@example creating_geometry
ring3 = GI.LinearRing(Point.(zip(x, y)))
````

Now create a `Polygon` from the `LineRing` with new the `CRS`

````@example creating_geometry
polygon3 = GI.Polygon([ring3]; crs=crs2)
````

Now plot on existing GeoAxis. NOTE: keyword argument `source` is used to specify the new `CRS`

````@example creating_geometry
plot!(ga,polygon3; color=:green, source=crs2)
fig
````

Great, we can make geometries and plot them on a map... now let's export the data to common geospatial data formats.

Typically, you'll also want to include attibutes with your geometries. Attibutes are simply data that is attibuted to each geometry. The easiest way to do this is to create a table with a `:geometry` column. Let's do this using [`DataFrames`](https://github.com/JuliaData/DataFrames.jl).

````@example creating_geometry
using DataFrames
import Shapefile
import GeoJSON
import GeoParquet


df = DataFrame(geometry=[polygon1, polygon2])
````

now lets add a couple of attributes to the geometries

````@example creating_geometry
df[!,:id] = ["a", "b"]
df[!, :name] = ["polygon 1", "polygon 2"]
````

now let's save as a `GeoJSON`

````@example creating_geometry
fn = "shapes.json"
GeoJSON.write(fn, df)
````

now let's save as a `Shapefile`

````@example creating_geometry
fn = "shapes.shp"
Shapefile.write(fn, df)
````

now let's save as a `GeoParquet`

````@example creating_geometry
fn = "shapes.parquet"
GeoParquet.write(fn, df, (:geometry,))
````

Finally, if there's no Julia-native package which saves in that format (`.gpkg`, `.gml`, etc), you can use [`GeoDataFrames`](https://github.com/evetion/GeoDataFrames.jl) to save your data.
This package uses the [GDAL](https://gdal.org/) library, which is the industry standard, via [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl).

````@example creating_geometry
import GeoDataFrames
fn = "shapes.gpkg"
GeoDataFrames.write(fn, df)
````

And there we go, you can now create mapped geometries from scratch, manipulate them, plot them on a map, and save
them in multiple geospatial data formats.
