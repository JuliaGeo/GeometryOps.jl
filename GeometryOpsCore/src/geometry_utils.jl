
# Trait-based _linearring to handle any GeoInterface-compatible geometry
_linearring(geom) = _linearring(GI.trait(geom), geom)

# If it's already a LinearRing, return as-is
_linearring(::GI.LinearRingTrait, geom) = geom

# If it's a LineString (e.g., from ArchGDAL), convert to LinearRing preserving CRS
_linearring(::GI.LineStringTrait, geom) =
    GI.LinearRing(GI.getpoint(geom); crs=GI.crs(geom))

# Concrete type specializations for GI wrappers (optimization)
_linearring(geom::GI.LineString) = GI.LinearRing(parent(geom); extent=geom.extent, crs=geom.crs)
_linearring(geom::GI.LinearRing) = geom
