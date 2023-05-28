
_is3d(geom) = _is3d(GI.trait(geom), geom)
_is3d(::GI.AbstractGeometryTrait, geom) = GI.is3d(geom)
_is3d(::GI.FeatureTrait, feature) = _is3d(GI.geometry(feature))
_is3d(::GI.FeatureCollectionTrait, fc) = _is3d(GI.getfeature(fc, 1))
_is3d(::Nothing, geom) = _is3d(first(geom)) # Otherwise step into an itererable
