# Geospatial Geometry

In this tutorial, we're going to: 
1. ...
2. ...
3. Geodesic paths...

First, we load some required packages: 

````julia
using Pkg
Pkg.add([])
````

````@example geospatial_geometry
# import ... 
````

Todo
Add packages
import
Change data example to be independent of creating_geometry.md


## [Plot geometries on a map using `GeoMakie` and coordinate reference system (`CRS`)](@id plot-geometry) 

In geospatial sciences we often have data in one [Coordinate Reference System (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) (`source`) and would like to display it in different (`destination`) `CRS`. `GeoMakie` allows us to do this by automatically projecting from `source` to `destination` CRS.

Here, our `source` CRS is common geographic (i.e. coordinates of latitude and longitude), [WGS84](https://epsg.io/4326).

````@example creating_geometry
source_crs1 = GFT.EPSG(4326)
````

Now let's pick a `destination` CRS for displaying our map. Here we'll pick [natearth2](https://proj.org/en/9.4/operations/projections/natearth2.html).

````@example creating_geometry
destination_crs = "+proj=natearth2"
````

Let's add land area for context. First, download and open the [Natural Earth](https://www.naturalearthdata.com) global land polygons at 110 m resolution.`GeoMakie` ships with this particular dataset, so we will access it from there.

````@example creating_geometry
land_path = GeoMakie.assetpath("ne_110m_land.geojson")
````

!!! note
    Natural Earth has lots of other datasets, and there is a Julia package that provides an interface to it called [NaturalEarth.jl](https://github.com/JuliaGeo/NaturalEarth.jl).

Read the land `MultiPolygon`s as a `GeoJSON.FeatureCollection`.

````@example creating_geometry
land_geo = GeoJSON.read(land_path)
````

We then need to create a figure with a `GeoAxis` that can handle the projection between `source` and `destination` CRS. For GeoMakie, `source` is the CRS of the input and `dest` is the CRS you want to visualize in.

````@example creating_geometry
fig = Figure(size=(1000, 500));
ga = GeoAxis(
    fig[1, 1];
    source = source_crs1,
    dest = destination_crs,
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

Now let's plot a `Polygon` like before, but this time with a CRS that differs from our `source` data

````@example creating_geometry
plot!(multipolygon; color = :green)
fig
````

But what if we want to plot geometries with a different `source` CRS on the same figure?

To show how to do this let's create a geometry with coordinates in UTM (Universal Transverse Mercator) zone 10N [EPSG:32610](https://epsg.io/32610).

````@example creating_geometry
source_crs2 = GFT.EPSG(32610)
````

Create a polygon (we're working in meters now, not latitude and longitude)
````@example creating_geometry
r = 1000000;
ϴ = 0:0.01:2pi;
x = r .* cos.(ϴ).^3 .+ 500000;
y = r .* sin.(ϴ) .^ 3 .+5000000;
DisplayAs.setcontext(y, :compact => true, :displaysize => (10, 40),) # hide
````

Now create a `LinearRing` from `Points`

````@example creating_geometry
ring3 = GI.LinearRing(Point.(zip(x, y)))
````

Now create a `Polygon` from the `LineRing` 

````@example creating_geometry
polygon3 = GI.Polygon([ring3])
````

Now plot on the existing GeoAxis. 

!!! note
    The keyword argument `source` is used to specify the source `CRS` of that particular plot, when plotting on an existing `GeoAxis`.

````@example creating_geometry
plot!(ga,polygon3; color=:red, source = source_crs2)
fig
````

## [Create geospatial geometries with embedded coordinate reference system information](@id geom-crs)
 
Great, we can make geometries and plot them on a map... now let's export the data to common geospatial data formats. To do this we now need to create geometries with embedded `CRS` information, making it a geospatial geometry. All that's needed is to include `; crs = crs` as a keyword argument when constructing the geometry.

Let's do this for a new `Polygon`
````@example creating_geometry
r = 3;
k = 7;
ϴ = 0:0.01:2pi;
x = r .* (k + 1) .* cos.(ϴ) .- r .* cos.((k + 1) .* ϴ);
y = r .* (k + 1) .* sin.(ϴ) .- r .* sin.((k + 1) .* ϴ);
ring4 = GI.LinearRing(Point.(zip(x, y)))
````

But this time when we create the `Polygon` we need to specify the `CRS` at the time of creation, making it a geospatial polygon

````@example creating_geometry
geopoly1 = GI.Polygon([ring4], crs = source_crs1)
````

!!! note
    It is good practice to only include CRS information with the highest-level geometry. Not doing so can bloat the memory footprint of the geometry. CRS information _can_ be included at the individual `Point` level but is discouraged.

And let's create second `Polygon` by shifting the first using CoordinateTransformations

````@example creating_geometry
xoffset = 20.;
yoffset = -25.;
f = CoordinateTransformations.Translation(xoffset, yoffset);
geopoly2 = GO.transform(f, geopoly1);
````

## [Creating a table with attributes and geometry](@id attributes)

Typically, you'll also want to include attributes with your geometries. Attributes are simply data that are attributed to each geometry. The easiest way to do this is to create a table with a `:geometry` column. Let's do this using [`DataFrames`](https://github.com/JuliaData/DataFrames.jl).

````@example creating_geometry
using DataFrames
df = DataFrame(geometry=[geopoly1, geopoly2])
````

Now let's add a couple of attributes to the geometries.  We do this using [DataFrames' `!` mutation syntax](https://dataframes.juliadata.org/stable/man/getting_started/#The-DataFrame-Type) that allows you to add a new column to an existing data frame.

````@example creating_geometry
df[!,:id] = ["a", "b"]
df[!, :name] = ["polygon 1", "polygon 2"]
df
````

## [Saving your geospatial data](@id save-geometry)

There are Julia packages for most commonly used geographic data formats.  Below, we show how to export that data to each of these.

We begin with [GeoJSON](https://github.com/JuliaGeo/GeoJSON.jl), which is a [JSON](https://en.wikipedia.org/wiki/JSON) format for geospatial feature collections.  It's human-readable and widely supported by most web-based and desktop geospatial libraries.

````@example creating_geometry
import GeoJSON
fn = "shapes.json"
GeoJSON.write(fn, df)
````

Now, let's save as a [`Shapefile`](https://github.com/JuliaGeo/Shapefile.jl).  Shapefiles are actually a set of files (usually 4) that hold geometry information, a CRS, and additional attribute information as a separate table.  When you give `Shapefile.write` a file name, it will write 4 files of the same name but with different extensions.

````@example creating_geometry
import Shapefile
fn = "shapes.shp"
Shapefile.write(fn, df)
````

Now, let's save as a [`GeoParquet`](https://github.com/JuliaGeo/GeoParquet.jl).  GeoParquet is a geospatial extension to the [Parquet](https://parquet.apache.org/) format, which is a high-performance data store.  It's great for storing large amounts of data in a single file.

````@example creating_geometry
import GeoParquet
fn = "shapes.parquet"
GeoParquet.write(fn, df, (:geometry,))
````

Finally, if there's no Julia-native package that can write data to your desired format (e.g. `.gpkg`, `.gml`, etc), you can use [`GeoDataFrames`](https://github.com/evetion/GeoDataFrames.jl). This package uses the [GDAL](https://gdal.org/) library under the hood which supports writing to nearly all geospatial formats.

````@example creating_geometry
import GeoDataFrames
fn = "shapes.gpkg"
GeoDataFrames.write(fn, df)
````

And there we go, you can now create mapped geometries from scratch, manipulate them, plot them on a map, and save them in multiple geospatial data formats.

## Geodesic paths

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