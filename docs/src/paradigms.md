# Paradigms

GeometryOps exposes functions like `apply` and `applyreduce`, as well as the `fix` and `prepare` APIs, that represent _paradigms_ of programming, by which we mean the ability to program in a certain way, and in so doing, fit neatly into the tools we've built without needing to re-implement the wheel.

Below, we'll describe some of the foundational paradigms of GeometryOps, and why you should care!

## `apply`


The `apply` function allows you to decompose a given collection of geometries down to a certain level, and then operate on it.  In general, its invocation is:

```julia
apply(f, trait::Trait, geom)
```

Functionally, it's similar to `map` in the way you apply it to geometries - except that you tell it at which level it should stop, by passing a `trait` to it.  

`apply` will start by decomposing the geometry, feature, featurecollection, iterable, or table that you pass to it, and stop when it encounters a geometry for which `GI.trait(geom) isa Trait`.  This encompasses unions of traits especially, but beware that any geometry which is not explicitly handled, and hits `GI.PointTrait`, will cause an error.

`apply` is unlike `map` in that it returns reconstructed geometries, instead of the raw output of the function.  If you want a purely map-like behaviour, like calculating the length of each linestring in your feature collection, then call  `GO.flatten(f, trait, geom)`, which will decompose each geometry to the given `trait` and apply `f` to it, returning the decomposition as a flattened vector.

### `applyreduce`

`applyreduce` is like the previous `map`-based approach that we mentioned, except that it `reduce`s the result of `f` by `op`.  Note that `applyreduce` does not guarantee associativity, so it's best to have `typeof(init) == returntype(op)`.

## `fix` and `prepare`

The `fix` and `prepare` paradigms are different from `apply`, though they are built on top of it.  They involve the use of structs as "actions", where a constructed object indicates an action that should be taken.  A trait like interface prescribes the level (polygon, linestring, point, etc) at which each action should be applied.

In general, the idea here is to be able to invoke several actions efficiently and simultaneously, for example when correcting invalid geometries, or instantiating a `Prepared` geometry with several preparations (sorted edge lists, rtrees, monotone chains, etc.)

