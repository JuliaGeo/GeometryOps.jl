#=

# Operations

!!! warning
    Operations are not yet implemented.  If you want to implement them then you may do so at your own risk - or file a PR!

Operations are callable structs, that contain the entire specification for what the algorithm will do.

Sometimes they may be underspecified and only materialized fully when you see the geometry, so you can extract
the best manifold for those geometries.

* Some conceptual thing that you do to a geometry
* Overloads on abstract type to decompose user input to have materialized algorithm, manifold, and geoms
* Run Operation{Alg{Manifold}}(trait, geom) at the lowest level
* Some indication on whether to use apply or applyreduce?  Or are we going too far here
    * if we do this, then we also need `operation_level` to return a geometry trait or traittarget

Operations may look like:

```julia
Arclength()(geoms)
Arclength(Geodesic())(geoms)
Arclength(Proj())(geoms)
Arclength(Proj(Geodesic(; ...)))(geoms)
Arclength(Ericsson())(geoms) # more precise, goes wonky if any points in a triangle are antipodal
Arclength(LHuilier())(geoms) # less precise, does not go wonky on antipodal points
```

Two argument operations, like polygon set operations, may look like:

```julia
Union(intersection_alg(manifold); exact, target)(geom1, geom2)
```

Here intersection_alg can be Foster, which we already have in GeometryOps, or GEOS
but if we ever implement e.g. RelateNG in GeometryOps, we can add that in.
=#

"""
    abstract type Operation{Alg <: Algorithm} end

Operations are callable structs, that contain the entire specification for what the algorithm will do.

Sometimes they may be underspecified and only materialized fully when you see the geometry, so you can extract
the best manifold for those geometries.
"""
abstract type Operation{Alg <: Algorithm} end

#= 
Here's an example of how this might work:

```julia
struct XPlusOneOperation{M <: Manifold} <: Operation{NoAlgorithm{M}}
    m::M # the manifold always needs to be stored, since it's not a singleton
    x::Int
end
```

```julia
struct Area{Alg<: Algorithm} <: Operation{Alg}
    alg::Alg
end

Area() = Area(NoAlgorithm(Planar()))

function (op::Area{Alg})(data; threaded = _False(), init = 0.0) where {Alg <: Algorithm}
    return GO.area(op.alg, data; threaded, init)
end
```

=#