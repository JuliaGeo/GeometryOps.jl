
_linearring(geom::GI.LineString) = GI.LinearRing(parent(geom); extent=geom.extent, crs=geom.crs)
_linearring(geom::GI.LinearRing) = geom
