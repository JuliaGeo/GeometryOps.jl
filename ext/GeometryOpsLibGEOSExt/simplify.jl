# Address potential ambiguities
GO._simplify(::GI.PointTrait, ::GO.TopologyPreserve, geom; kw...) = geom
GO._simplify(::GI.MultiPointTrait, ::GO.TopologyPreserve, geom; kw...) = geom

function GO._simplify(::GI.AbstractGeometryTrait, alg::GO.GEOS, geom)
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

