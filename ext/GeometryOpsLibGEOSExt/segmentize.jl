# # Segmentize
import GeometryOps: segmentize, apply

#=
This file implements the LibGEOS segmentization method for GeometryOps.
=#

function _segmentize_geos(geom::LG.AbstractGeometry, max_distance)
    context = LG.get_context(geom)
    result = LG.GEOSDensify_r(context, geom, max_distance)
    if result == C_NULL
        error("LibGEOS: Error in GEOSDensify")
    end
    return LG.geomFromGEOS(result, context)
end

_segmentize_geos(geom, max_distance) = _segmentize_geos(GI.convert(LG, geom), max_distance)

function GO.segmentize(alg::GEOS{(:max_distance,)}, geom; threaded::Union{Bool, GO.BoolsAsTypes} = _False())
    return GO.apply(
        Base.Fix2(_segmentize_geos, alg.params.max_distance), 
        GO.TraitTarget(GI.GeometryCollectionTrait(), GI.MultiPolygonTrait(), GI.PolygonTrait(), GI.MultiLineStringTrait(), GI.LineStringTrait(), GI.LinearRingTrait(), GI.MultiPointTrait(), GI.PointTrait()),
        geom; 
        threaded
    )
end