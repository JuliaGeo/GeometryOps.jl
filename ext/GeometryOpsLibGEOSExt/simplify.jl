# Address potential ambiguities
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

