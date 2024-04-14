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

# 2 behaviours:
# - enforce: enforce the presence of a kwargs
# - fetch: fetch the value of a kwargs, or return a default value
function GO.segmentize(alg::GEOS, geom; threaded::Union{Bool, GO.BoolsAsTypes} = _False())
    max_distance = enforce(alg, :max_distance, GO.segmentize)
    return GO.apply(
        Base.Fix2(_segmentize_geos, max_distance), 
        GO.TraitTarget(GI.GeometryCollectionTrait(), GI.MultiPolygonTrait(), GI.PolygonTrait(), GI.MultiLineStringTrait(), GI.LineStringTrait(), GI.LinearRingTrait(), GI.MultiPointTrait(), GI.PointTrait()),
        geom; 
        threaded
    )
end