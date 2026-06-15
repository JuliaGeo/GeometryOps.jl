# Geospatial Geometry

In this tutorial, we're going to: 
1. [Plot geometries on a map using `GeoMakie` and coordinate reference system (`CRS`)](@ref plot-geometry)
2. [Create geospatial geometries with embedded coordinate reference system information](@ref geom-crs)
3. [Assign attributes to geospatial geometries](@ref attributes)
4. [Save geospatial geometries to common geospatial file formats](@ref save-geometry)
5. [Introduce Geodesic Paths](@ref geodesic-paths)

Install the packages used in this tutorial:

````julia
using Pkg
Pkg.add(["GeoInterface", "GeometryOps", "GeoFormatTypes", 
        "GeoJSON", "GeoParquet", "GeoDataFrames",
        "CoordinateTransformations", "Proj", "DataFrames", 
        "CairoMakie", "GeoMakie", "Shapefile"])
````

````@example geospatial_geometry
# Geospatial packages from Julia
import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
using GeoParquet
using GeoDataFrames
# Packages for coordinate transformation and projection
import CoordinateTransformations
import Proj
# Plotting
using CairoMakie
using GeoMakie
using DisplayAs # hide
Makie.set_theme!(Makie.MAKIE_DEFAULT_THEME) # hide
# Loading data
using GeoJSON 
using DataFrames
import Shapefile
````


## [Plot geometries on a map using `GeoMakie` and coordinate reference system (`CRS`)](@id plot-geometry) 

In geospatial sciences we often have data in one [Coordinate Reference System (CRS)](https://en.wikipedia.org/wiki/Spatial_reference_system) (`source`) and would like to display it in different (`destination`) `CRS`. `GeoMakie` allows us to do this by automatically projecting from `source` to `destination` CRS.

Here, our `source` CRS is common geographic (i.e. coordinates of latitude and longitude), [WGS84](https://epsg.io/4326).

````@example geospatial_geometry
source_crs1 = GFT.EPSG(4326)
````

Now let's pick a `destination` CRS for displaying our map. Here we'll pick [natearth2](https://proj.org/en/9.4/operations/projections/natearth2.html).

````@example geospatial_geometry
destination_crs = "+proj=natearth2"
````

Let's add land area for context. First, download and open the [Natural Earth](https://www.naturalearthdata.com) global land polygons at 110 m resolution.  `GeoMakie` ships with this particular dataset, so we will access it from there.
Note that this will be a path on your local machine, so you could easily point it to any other `.geojson` file you have.

````@example geospatial_geometry
land_path = GeoMakie.assetpath("ne_110m_land.geojson")
````

!!! note
    Natural Earth has lots of other datasets, and there is a Julia package that provides an interface to it called [NaturalEarth.jl](https://github.com/JuliaGeo/NaturalEarth.jl).

Read this dataset in as a `GeoJSON.FeatureCollection`.

````@example geospatial_geometry
land_geo = GeoJSON.read(land_path)
````

We then need to create a figure with a `GeoAxis` that can handle the projection between `source` and `destination` CRS. For [`GeoMakie`](https://geo.makie.org/stable/), `source` is the CRS of the input and `dest` is the CRS you want to visualize in.

````@example geospatial_geometry
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

````@example geospatial_geometry
poly!(ga, land_geo; color=:black)
fig
````

Now let's plot a `Polygon` like before, but this time with a CRS that differs from our `source` data

````@example geospatial_geometry
function ring(radius)
    ϴ = 0:0.01:2π
	points = @. GI.Point(50 + radius * cos(ϴ), 50 + radius * sin(ϴ))
	return GI.LinearRing(points)
end

function spiro(a = 22, b = 2, k = 11)
    ϴ = 0:0.01:2π
    points = @. GI.Point(
        50 + a*cos(ϴ) - b*cos(k*ϴ),
        50 + a*sin(ϴ) - b*sin(k*ϴ)
    )
    return GI.LinearRing(points)
end


multipolygon = GI.MultiPolygon([
    GI.Polygon([spiro(), ring(8)]),
    GI.Polygon([ring(4)])
])

plot!(multipolygon; color = :green)
fig
````

But what if we want to plot geometries with a different `source` CRS on the same figure?

To show how to do this let's create a geometry with coordinates in UTM (Universal Transverse Mercator) zone 10N [EPSG:32610](https://epsg.io/32610).

````@example geospatial_geometry
source_crs2 = GFT.EPSG(32610)
````

Create a polygon (we're working in meters now, not latitude and longitude)
````@example geospatial_geometry
r = 1000000;
ϴ = 0:0.01:2pi;
x = r .* cos.(ϴ).^3 .+ 500000;
y = r .* sin.(ϴ) .^ 3 .+5000000;
DisplayAs.setcontext(y, :compact => true, :displaysize => (10, 40),) # hide
````

Now create a `LinearRing` from `Points`

````@example geospatial_geometry
ring3 = GI.LinearRing(Point.(zip(x, y)))
````

Now create a `Polygon` from the `LineRing` 

````@example geospatial_geometry
polygon3 = GI.Polygon([ring3])
````

Now plot on the existing GeoAxis. 

!!! note
    The keyword argument `source` is used to specify the source `CRS` of that particular plot, when plotting on an existing `GeoAxis`.

````@example geospatial_geometry
plot!(ga,polygon3; color=:red, source = source_crs2)
fig
````

## [Create geospatial geometries with embedded coordinate reference system information](@id geom-crs)
 
Great, we can make geometries and plot them on a map... now let's export the data to common geospatial data formats. To do this we now need to create geometries with embedded `CRS` information, making it a geospatial geometry. All that's needed is to include `; crs = crs` as a keyword argument when constructing the geometry.

Let's do this for a new `Polygon`
````@example geospatial_geometry
r = 3;
k = 7;
ϴ = 0:0.01:2pi;
x = r .* (k + 1) .* cos.(ϴ) .- r .* cos.((k + 1) .* ϴ);
y = r .* (k + 1) .* sin.(ϴ) .- r .* sin.((k + 1) .* ϴ);
ring4 = GI.LinearRing(Point.(zip(x, y)))
````

But this time when we create the `Polygon` we need to specify the `CRS` at the time of creation, making it a geospatial polygon

````@example geospatial_geometry
geopoly1 = GI.Polygon([ring4], crs = source_crs1)
````

!!! note
    It is good practice to only include CRS information with the highest-level geometry. Not doing so can bloat the memory footprint of the geometry. CRS information _can_ be included at the individual `Point` level but is discouraged.

And let's create a second `Polygon` by shifting the first using CoordinateTransformations

````@example geospatial_geometry
xoffset = 20.
yoffset = -25.
f = CoordinateTransformations.Translation(xoffset, yoffset)
geopoly2 = GO.transform(f, geopoly1)
````

## [Creating a table with attributes and geometry](@id attributes)

Typically, you'll also want to include attributes with your geometries. Attributes are simply data that are attributed to each geometry. The easiest way to do this is to create a table with a `:geometry` column. Let's do this using [`DataFrames`](https://github.com/JuliaData/DataFrames.jl).

````@example geospatial_geometry
df = DataFrame(geometry=[geopoly1, geopoly2])
````

Now let's add a couple of attributes to the geometries by adding new columns to our existing data frame.

````@example geospatial_geometry
df.id = ["a", "b"]
df.name = ["polygon 1", "polygon 2"]
df
````

## [Saving your geospatial data](@id save-geometry)

There are Julia packages for most commonly used geographic data formats.  Below, we show how to export that data to each of these.

::: tabs

== GeoDataFrames

In general, [`GeoDataFrames`](https://github.com/evetion/GeoDataFrames.jl) is recomended as the default way to write data to your desired format. This package uses the [GDAL](https://gdal.org/) library under the hood which supports writing to nearly all geospatial formats.

Writing to `gpkg`:

````@example geospatial_geometry
GeoDataFrames.write("shapes.gpkg", df)
````

Writing to `GeoJSON`:

````@example geospatial_geometry
GeoDataFrames.write("file.geojson", df)
````

View the [`GeoDataFrames`](https://github.com/evetion/GeoDataFrames.jl) documentation for all recognized file extensions.


== GeoJSON

Next, let's save as a [GeoJSON](https://github.com/JuliaGeo/GeoJSON.jl), which is a [JSON](https://en.wikipedia.org/wiki/JSON) format for geospatial feature collections.  It's human-readable and widely supported by most web-based and desktop geospatial libraries.

````@example geospatial_geometry
GeoJSON.write("shapes.json", df)
````

== Shapefile

Now, let's save as a [`Shapefile`](https://github.com/JuliaGeo/Shapefile.jl).  Shapefiles are actually a set of files (usually 4) that hold geometry information, a CRS, and additional attribute information as a separate table.  When you give `Shapefile.write` a file name, it will write 4 files of the same name but with different extensions.

````@example geospatial_geometry
Shapefile.write("shapes.shp", df)
````

== GeoParquet

Now, let's save as a [`GeoParquet`](https://github.com/JuliaGeo/GeoParquet.jl).  GeoParquet is a geospatial extension to the [Parquet](https://parquet.apache.org/) format, which is a high-performance data store.  It's great for storing large amounts of data in a single file.

````@example geospatial_geometry
GeoParquet.write("shapes.parquet", df)
````

:::

## [Geodesic paths](@id geodesic-paths)

Geodesic paths are paths computed on an ellipsoid, as opposed to a plane. The geodesic is the shortest path between two points measured along the Earth's curved surface. Because the surface is curved, that shortest path appears as a curve (not a straight line) when drawn on a flat map. 

Here, we use the `segmentize` function to add vertices along the geodesic between two points, so the line follows Earth's curved surface instead of a straight line. 

````@example geospatial_geometry
# Two points in (longitude, latitude) order: Houston (IAH) and Amsterdam (AMS).
# Geodesic methods assume lon/lat input.
IAH = (-95.358421, 29.749907)
AMS = (4.897070, 52.377956)

# Draw coastlines for geographic context
fig, ga, _cp = lines(GeoMakie.coastlines(); axis = (; type = GeoAxis))

# Create our line along the Earth, accounting for curvature
lines!(ga, GO.segmentize(GO.Geodesic(), GI.LineString([IAH, AMS]); max_distance = 100_000))
fig
````

And there we go, you can now create mapped geometries from scratch, manipulate them, plot them on a map, and save them in multiple geospatial data formats, as well as create geodesic paths. 
