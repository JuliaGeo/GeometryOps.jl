# Peculiarities

## What does `apply` return and why?

`apply` returns the target geometries returned by `f`, whatever type/package they are from, but geometries, features or feature collections that wrapped the target are replaced with GeoInterace.jl wrappers with matching `GeoInterface.trait` to the originals. All non-geointerface iterables become `Array`s. Tables.jl compatible tables are converted either back to the original type if a `Tables.materializer` is defined, and if not then returned as generic `NamedTuple` column tables (i.e., a NamedTuple of vectors).

 It is recommended for consistency that `f` returns GeoInterface geometries unless there is a performance/conversion overhead to doing that. 

## Why do you want me to provide a `target` in set operations?

In polygon set operations like `intersection`, `difference`, and `union`, many different geometry types may be obtained - depending on the relationship between the polygons.  For example, when performing an union on two nonintersecting polygons, one would technically have two disjoint polygons as an output.

We use the `target` keyword to allow the user to control which kinds of geometry they want back.  For example, setting `target` to `PolygonTrait` will cause a vector of polygons to be returned (this is the only currently supported behaviour).  In future, we may implement `MultiPolygonTrait` or `GeometryCollectionTrait` targets which will return a single geometry, as LibGEOS and ArchGDAL do.

This also allows for a lot more type stability - when you ask for polygons, we won't return a geometrycollection with line segments.  Especially in simulation workflows, this is excellent for simplified data processing.

## `True` and `False` (or `BoolsAsTypes`)

!!! warning
    These are internals and explicitly *not* public API,
    meaning they may change at any time!

When dispatch can be controlled by the value of a boolean variable, this introduces type instability.  Instead of introducing type instability, we chose to encode our boolean decision variables, like `threaded` and `calc_extent` in `apply`, as types.  This allows the compiler to reason about what will happen, and call the correct compiled method, in a stable way without worrying about 