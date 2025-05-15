#=
# GEOS simplify

This file implements [`GO.simplify`](@ref) using the [`GEOS`](@ref) algorithm.

You can pass the keyword `method` to GEOS, to choose between the `:TopologyPreserve` and `:DouglasPeucker` methods.

The only other parameter required is `tol` which defines the simplification tolerance.

## Example

```@example geos-simplify
import GeometryOps as GO, GeoInterface as GI
import LibGEOS # activate the GEOS algorithm
using GADM # get data
using CairoMakie # to plot

france = GI.geometry(GI.getfeature(GADM.get("FRA"), 1))

simplified = GO.simplify(GO.GEOS(; tol = 1), france)

f, a, p = lines(france; label = "Original")
lines!(a, simplified; label = "Simplified")
axislegend(a)
f
```

You can also choose the TopologyPreserve method:

```@example geos-simplify
simplified = GO.simplify(GO.GEOS(; tol = 1, method = :TopologyPreserve), france)

f, a, p = lines(france; label = "Original")
lines!(a, simplified; label = "Topology preserving simplify")
axislegend(a)
f
```

=#

function GO.simplify(alg::GO.GEOS, data; kwargs...)
    method = get(alg, :method, :TopologyPreserve)
    tol = enforce(alg, :tol, GO.simplify)

    functor = if method == :TopologyPreserve
        (trait, geom) -> LG.topologyPreserveSimplify(GI.convert(LG.geointerface_geomtype(trait), trait, geom), tol)
    elseif method == :DouglasPeucker
        (trait, geom) -> LG.simplify(GI.convert(LG.geointerface_geomtype(trait), trait, geom), tol)
    else
        error("Invalid method passed to `GO.simplify(GEOS(...), ...)`: $method. \nPlease use :TopologyPreserve or :DouglasPeucker")
    end

    return GO.apply(WithTrait(functor), GO.TraitTarget{GI.AbstractGeometryTrait}(), data; kwargs...)
end
