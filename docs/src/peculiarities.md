# Peculiarities

## `_True` and `_False` (or `BoolsAsTypes`)

When dispatch can be controlled by the value of a boolean variable, this introduces type instability.  Instead of introducing type instability, we chose to encode our boolean decision variables, like `threaded` and `calc_extent` in `apply`, as types.  This allows the compiler to reason about what will happen, and call the correct compiled method, in a stable way without worrying about 

## What does `apply` return and why?

`apply` returns the target geometries returned by `f`, whatever type/package they are from, but geometries, features or feature collections that wrapped the target are replaced with GeoInterace.jl wrappers with matching `GeoInterface.trait` to the originals. All non-geointerface iterables become `Array`s. Tables.jl compatible tables are converted either back to the original type if a `Tables.materializer` is defined, and if not then returned as generic `NamedTuple` column tables (i.e., a NamedTuple of vectors).

 It is recommended for consistency that `f` returns GeoInterface geometries unless there is a performance/conversion overhead to doing that. 

## Why do you want me to provide a `target` in set operations?

@skygering

Mainly type stability reasons.
