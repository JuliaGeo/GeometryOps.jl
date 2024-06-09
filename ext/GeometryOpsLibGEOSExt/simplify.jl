# Address potential ambiguities
GO._simplify(::GI.PointTrait, ::GO.TopologyPreserve, geom; kw...) = geom
GO._simplify(::GI.MultiPointTrait, ::GO.TopologyPreserve, geom; kw...) = geom

function GO._simplify(::GI.AbstractGeometryTrait, alg::GO.TopologyPreserve, geom)
    return LG.topologyPreserveSimplify(GI.convert(LG, geom), alg.tol)
end

