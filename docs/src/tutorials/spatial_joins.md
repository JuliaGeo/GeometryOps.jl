# Spatial joins

Spatial joins are joins which are based not on equality, but on some predicate ``p(x, y)``, which takes two geometries, and returns a value of either `true` or `false`.  For geometries, the [`DE-9IM`](https://en.wikipedia.org/wiki/DE-9IM) spatial relationship model is used to determine the spatial relationship between two geometries.  

In this tutorial, we will show how to perform a spatial join on first a toy dataset and then two Natural Earth datasets, to show how this can be used in the real world.

In order to perform the spatial join, we use [FlexiJoins.jl](https://github.com/JuliaAPlavin/FlexiJoins.jl) to perform the join, specifically using its `by_pred` joining method.  This allows the user to specify a predicate in the following manner:
```julia
[inner/left/outer/...]join((table1, table1),
    by_pred(:table1_column, predicate_function, :table2_column)
)
```

We have enabled the use of all of GeometryOps' boolean comparisons here.  These are:

```julia
GO.contains, GO.within, GO.intersects, GO.touches, GO.crosses, GO.disjoint, GO.overlaps, GO.covers, GO.coveredby, GO.equals
```

## Simple example

This example demonstrates how to perform a spatial join between two datasets: a set of polygons and a set of randomly generated points. 

The polygons are represented as a DataFrame with geometries and colors, while the points are stored in a separate DataFrame. 

The spatial join is performed using the `contains` predicate from GeometryOps, which checks if each point is contained within any of the polygons. The resulting joined DataFrame is then used to plot the points, colored according to the containing polygon.

First, we generate our data.  We create two triangle polygons which, together, span the rectangle (0, 0, 1, 1), and a set of points which are randomly distributed within this rectangle.

```@example spatialjoins
import GeoInterface as GI, GeometryOps as GO
using FlexiJoins, DataFrames

using CairoMakie, GeoInterfaceMakie

pl = GI.Polygon([GI.LinearRing([(0, 0), (1, 0), (1, 1), (0, 0)])])
pu = GI.Polygon([GI.LinearRing([(0, 0), (0, 1), (1, 1), (0, 0)])])
poly_df = DataFrame(geometry = [pl, pu], color = [:red, :blue])
f, a, p = Makie.with_theme(Attributes(; Axis = (; aspect = DataAspect()))) do # hide
f, a, p = poly(poly_df.geometry; color = tuple.(poly_df.color, 0.3))
end # hide
```

Here, the upper polygon is blue, and the lower polygon is red.  Keep this in mind!

Now, we generate the points.

```@example spatialjoins
points = tuple.(rand(100), rand(100))
points_df = DataFrame(geometry = points)
scatter!(points_df.geometry)
f
```

You can see that they are evenly distributed around the box.  But how do we know which points are in which polygons?

The answer here is to perform a spatial join.

Now, we can perform the "spatial join" using FlexiJoins.  We are performing an outer join here

```@example spatialjoins
joined_df = FlexiJoins.innerjoin(
    (poly_df, points_df), 
    by_pred(:geometry, GO.contains, :geometry)
)
```

```@example spatialjoins
scatter(joined_df.geometry_1; color = joined_df.color)
```

Here, you can see that the colors were assigned appropriately to the scattered points!

## Real-world example

Suppose I have a list of polygons representing administrative regions (or mining sites, or what have you), and I have a list of polygons for each country.  I want to find the country each region is in.

