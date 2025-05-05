# # GEOS segmentize
# This file implements [`segmentize`](@ref) using the [`GEOS`](@ref) algorithm.
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

function _wrap_and_segmentize_geos(geom, max_distance)
    _wrap(_segmentize_geos(geom, max_distance); crs = GI.crs(geom), calc_extent = false)
end

# 2 behaviours:
# - enforce: enforce the presence of a kwargs
# - fetch: fetch the value of a kwargs, or return a default value
@inline function GO.segmentize(alg::GEOS, geom; threaded::Union{Bool, GO.BoolsAsTypes} = False())
    max_distance = enforce(alg, :max_distance, GO.segmentize)
    return GO.apply(
        Base.Fix2(_wrap_and_segmentize_geos, max_distance), 
        # TODO: should this just be a target on GI.AbstractGeometryTrait()?
        # But Geos doesn't support eg RectangleTrait
        # Maybe we need an abstract trait `GI.AbstractWKBGeomTrait`?
        GO.TraitTarget(GI.GeometryCollectionTrait(), GI.MultiPolygonTrait(), GI.PolygonTrait(), GI.MultiLineStringTrait(), GI.LineStringTrait(), GI.LinearRingTrait(), GI.MultiPointTrait(), GI.PointTrait()),
        geom; 
        threaded
    )
end