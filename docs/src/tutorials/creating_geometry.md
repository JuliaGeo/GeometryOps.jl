# Creating geometry

````@example creating_geometry
import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import CoordinateTransformations
import Proj
using CairoMakie
````

The first thing we need is to decide what Coodinate reference system (crs) that we'll be
working in. Here we start with the most common geographic crs [WGS84](https://epsg.io/4326)

````@example creating_geometry
crs = GFT.EPSG(4326)
````

Let's start by making a single point with crs info included

````@example creating_geometry
point = GI.Point(0, 0; crs)
````

now let's plot our point

````@example creating_geometry
fig, ax, plt = plot(point)
````

now let's create a set of points, and have a bit more fun with plotting

````@example creating_geometry
x = [-5, 0, 5, 0]
y = [0, -5, 0, 5]
points = GI.Point.(zip(x,y); crs)
plot!(ax, points; marker = '✈', markersize = 30)
fig
````

points can be combined into a single MultiPoint geometry

````@example creating_geometry
x = [-5, -5, 5, 5]
y = [-5, 5, 5, -5]
multipoint = GeoInterface.MultiPoint(GI.Point.(zip(x, y); crs))
plot!(ax, multipoint, marker = '☁', markersize = 30)
display(fig)
````

let's create a line between two points

````@example creating_geometry
p1 = GI.Point.(-5, 0; crs)
p2 = GI.Point.(5, 0; crs)
line = GI.LineString([p1,p2]; crs)
plot!(ax, line)
display(fig)
````

now lets create a line connecting multiple points (i.e. a LineString)
this time getting a bit more fancy with point creation

````@example creating_geometry
r = 2;
k = 10;
ϴ = 0:0.01:2pi
x = r .* (k + 1) .* cos.(ϴ) .- r .* cos.((k + 1) .* ϴ)
y = r .* (k + 1) .* sin.(ϴ) .- r .* sin.((k + 1) .* ϴ)
lines = GI.LineString(GI.Point.(zip(x,y)); crs)
plot!(ax, lines; linewidth = 3)
display(fig)
````

We can also create a single LinearRing Trait, the building blocks of a polygon
A LiearRing is simply a LineString with the same begin and endpoint.
i.e. an arbitraty closed shape composed of point pairs

a LinearRing is composed of series of points listed in clockwise order (i.e. winding order)
I always think of a polygons as filled to the right of the lines as one progresses
from point n to point n+1

````@example creating_geometry
ring1 = GI.LinearRing(GI.getpoint(lines))
````

now lets make the LineRing into a Polygon

````@example creating_geometry
polygon1 = GI.Polygon([ring1]; crs)
````

now we can use GeometryOperations and CoordinateTransformations to shift polygon1
up to avoid plotting over our earlier results

````@example creating_geometry
xoffset = 0.
yoffset = 50.
f = CoordinateTransformations.Translation(xoffset, yoffset)
polygon1 = GO.transform(f, polygon1)
plot!(polygon1)
display(fig)
````

Polygons can contain "holes". The first LineRing in a polygon, all subsequent LineRings
are treated as holes in the leadind LineRing

NOTE: some packages consider the secondary LineRings holes, others look at the winding
order, where the polygons if filled inward if it has a clockwise winding order and outward
if it has a counterclockwise winding order... hopfully these are details that you'll never
have to deal with but are good to know

````@example creating_geometry
hole = GI.LinearRing(GI.getpoint(multipoint))
polygon1 = GI.Polygon([ring1, hole]; crs)
````

shift multiepolygon to the righ to avoid plotting over our earlier results

````@example creating_geometry
xoffset = 50.
yoffset = 0.
f = CoordinateTransformations.Translation(xoffset, yoffset)
polygon1 = GO.transform(f, polygon1)
plot!(polygon1)
display(fig)
````

Polygons can also be grouped together as a MultiPolygon

````@example creating_geometry
r = 5
x = cos.(reverse(ϴ)) .* r .+ xoffset
y = sin.(reverse(ϴ)) .* r .+ yoffset
ring2 =  GI.LinearRing(GI.Point.(zip(x,y)))
polygon2 = GI.Polygon([ring2])
multipolygon = GI.MultiPolygon([polygon1, polygon2]; crs)
````

shift multiepolygon up to avoid plotting over our earlier results

````@example creating_geometry
xoffset = 0.
yoffset = 50.
f = CoordinateTransformations.Translation(xoffset, yoffset)
multipolygon = GO.transform(f, multipolygon)
plot!(multipolygon)
display(fig)
````

Great now we can make Points, MultiPoints, Lines, LineStrings, Polygons (w holes), and
MultiPolygons.

But where the crs informatio come in? To show this we need to use GeoMakie that can
interpret the crs infomation that we've included with our geometries

add additional packages

````@example creating_geometry
using GeoMakie
using GeoMakie: GeoJSON
using Downloads
````

Now specify the source and destination projections for our map. Rememebr that the very
first thing we did was set our source coordinate system

````@example creating_geometry
source = crs;
dest = "+proj=natearth2" #see [https://proj.org/en/9.4/operations/projections/natearth2.html]
````

download Natural Earth continental outlines (https://www.naturalearthdata.com/)

````@example creating_geometry
land = GeoMakie.assetpath("ne_110m_land.geojson")
land_geo = GeoJSON.read(read(land, String))
````

create a figure with a GeoAxis that can handle the projections

````@example creating_geometry
fig = Figure(size=(1000, 500));
ga = GeoAxis(
    fig[1, 1];
    source=crs,
    dest=dest,
    xticklabelsvisible = false,
    yticklabelsvisible = false,
);
nothing #hide
````

plot land for context

````@example creating_geometry
poly!(ga, land_geo, color=:black)
display(fig)
````

now let's make a polygon like before

````@example creating_geometry
plot!(multipolygon; color = :green)
display(fig)
````

Great, we can make geometries and plot them on a map... not let's export the data to
common geospatial data formats

Typically you'll also want to include attibutes with your geometries. The easiest way to
do that is to create a table with a `:geometry` column. Let's do this using DataFrames

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

now let's save as a GeoJSON

````@example creating_geometry
fn = "shapes.json"
GeoJSON.write(fn, df)
````

now let's save as a Shapefile

````@example creating_geometry
fn = "shapes.shp"
Shapefile.write(fn, df)
````

now let's save as a GeoParquet

````@example creating_geometry
fn = "shapes.parquet"
GeoParquet.write(fn, df, (:geometry,))
````

and there we go, you can now create mapped geometries from scratch, plot on a map and save
in multiple geospatial data formats


