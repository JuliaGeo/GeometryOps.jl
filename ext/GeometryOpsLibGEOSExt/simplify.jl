# Address potential ambiguities
GO._simplify(::GI.PointTrait, ::GO.TopologyPreserve, geom; kw...) = geom
GO._simplify(::GI.MultiPointTrait, ::GO.TopologyPreserve, geom; kw...) = geom

function GO._simplify(::GI.AbstractGeometryTrait, alg::GO.TopologyPreserve, geom)
    return LIBGEOS.topologyPreserveSimplify(GI.convert(LibGEOS, geom), alg.tol)
end

