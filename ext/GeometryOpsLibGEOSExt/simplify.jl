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
GO._simplify(::GI.PointTrait, ::GO.GEOS, geom; kw...) = geom
GO._simplify(::GI.MultiPointTrait, ::GO.GEOS, geom; kw...) = geom

function GO._simplify(::GI.AbstractGeometryTrait, alg::GO.GEOS, geom; kwargs...)
    method = get(alg, :method, :TopologyPreserve)
    @assert haskey(alg.params, :tol) """
        The `:tol` parameter is required for the GEOS algorithm in `simplify`, 
        but it was not provided.  
        
        Provide it by passing `GEOS(; tol = ...,) as the algorithm.
        """
    tol = alg.params.tol
    if method == :TopologyPreserve
        return LG.topologyPreserveSimplify(GI.convert(LG, geom), tol)
    elseif method == :DouglasPeucker
        return LG.simplify(GI.convert(LG, geom), tol)
    else
        error("Invalid method passed to `GO.simplify(GEOS(...), ...)`: $method. Please use :TopologyPreserve or :DouglasPeucker")
    end
end

function GO._simplify(trait::GI.AbstractCurveTrait, alg::GO.GEOS, geom; kw...)
    Base.invoke(
        GO._simplify,
        Tuple{GI.AbstractGeometryTrait, GO.GEOS, typeof(geom)},
        trait, alg, geom; 
        kw...
    )
end

