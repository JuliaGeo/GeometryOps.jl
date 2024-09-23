# Spatial joins

Spatial joins are [table joins](https://www.geeksforgeeks.org/sql-join-set-1-inner-left-right-and-full-joins/) which are based not on equality, but on some predicate ``p(x, y)``, which takes two geometries, and returns a value of either `true` or `false`.  For geometries, the [`DE-9IM`](https://en.wikipedia.org/wiki/DE-9IM) spatial relationship model is used to determine the spatial relationship between two geometries.  

Spatial joins can be done between any geometry types (from geometrycollections to points), just as geometrical predicates can be evaluated on any geometries.

In this tutorial, we will show how to perform a spatial join on first a toy dataset and then two Natural Earth datasets, to show how this can be used in the real world.

In order to perform the spatial join, we use **[FlexiJoins.jl](https://github.com/JuliaAPlavin/FlexiJoins.jl)** to perform the join, specifically using its `by_pred` joining method.  This allows the user to specify a predicate in the following manner, for any kind of table join operation:
```julia
using FlexiJoins
innerjoin((table1, table1),
    by_pred(:table1_column, predicate_function, :table2_column) # & add other conditions here
)
leftjoin((table1, table1),
    by_pred(:table1_column, predicate_function, :table2_column) # & add other conditions here
)
rightjoin((table1, table1),
    by_pred(:table1_column, predicate_function, :table2_column) # & add other conditions here
)
outerjoin((table1, table1),
    by_pred(:table1_column, predicate_function, :table2_column) # & add other conditions here
)
```

We have enabled the use of all of GeometryOps' boolean comparisons here.  These are:

```julia
GO.contains, GO.within, GO.intersects, GO.touches, GO.crosses, GO.disjoint, GO.overlaps, GO.covers, GO.coveredby, GO.equals
```

!!! tip
    Always place the dataframe with more complex geometries second, as that is the one which will be sorted into a tree.

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
points = tuple.(rand(1000), rand(1000))
points_df = DataFrame(geometry = points)
scatter!(points_df.geometry)
f
```

You can see that they are evenly distributed around the box.  But how do we know which points are in which polygons?

We have to join the two dataframes based on which polygon (if any) each point lies within.

Now, we can perform the "spatial join" using FlexiJoins.  We are performing an outer join here

```@example spatialjoins
@time joined_df = FlexiJoins.innerjoin(
    (points_df, poly_df), 
    by_pred(:geometry, GO.within, :geometry)
)
```

```@example spatialjoins
scatter!(a, joined_df.geometry; color = joined_df.color)
f
```

Here, you can see that the colors were assigned appropriately to the scattered points!

## Real-world example

Suppose I have a list of polygons representing administrative regions (or mining sites, or what have you), and I have a list of polygons for each country.  I want to find the country each region is in.

```julia real
import GeoInterface as GI, GeometryOps as GO
using FlexiJoins, DataFrames, GADM # GADM gives us country and sublevel geometry

using CairoMakie, GeoInterfaceMakie

country_df = GADM.get.(["JPN", "USA", "IND", "DEU", "FRA"]) |> DataFrame
country_df.geometry = GI.GeometryCollection.(GO.tuples.(country_df.geom))

state_doublets = [
    ("USA", "New York"),
    ("USA", "California"),
    ("IND", "Karnataka"),
    ("DEU", "Berlin"),
    ("FRA", "Grand Est"),
    ("JPN", "Tokyo"),
]

state_full_df = (x -> GADM.get(x...)).(state_doublets) |> DataFrame
state_full_df.geom = GO.tuples.(only.(state_full_df.geom))
state_compact_df = state_full_df[:, [:geom, :NAME_1]]
```

```julia real
innerjoin((state_compact_df, country_df), by_pred(:geom, GO.within, :geometry))
innerjoin((state_compact_df,  view(country_df, 1:1, :)), by_pred(:geom, GO.within, :geometry))
```

!!! warning
    This is how you would do this, but it doesn't work yet, since the GeometryOps predicates are quite slow on large polygons.  If you try this, the code will continue to run for a very, very long time (it took 12 hours on my laptop, but with minimal CPU usage).

## Enabling custom predicates

In case you want to use a custom predicate, you only need to define a method to tell FlexiJoins how to use it.

For example, let's suppose you wanted to perform a spatial join on geometries which are some distance away from each other:

```julia
my_predicate_function = <(5) ∘ abs ∘ GO.distance
```

You would need to define `FlexiJoins.supports_mode` on your predicate:

```julia{3}
FlexiJoins.supports_mode(
    ::FlexiJoins.Mode.NestedLoopFast, 
    ::FlexiJoins.ByPred{typeof(my_predicate_function)}, 
    datas
) = true
```

This will enable FlexiJoins to support your custom function, when it's passed to `by_pred(:geometry, my_predicate_function, :geometry)`.
