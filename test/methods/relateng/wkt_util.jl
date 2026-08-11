# Shared WKT-loading helper for the relateng test files.
#
# Mirrors `jts_wkt_to_geom` in test/external/jts/jts_testset_reader.jl (which
# keeps its own copy for the XML harness): WellKnownGeometry → GO.tuples for
# plain WKT, with a LibGEOS fallback for WKT that WellKnownGeometry or GI
# wrapper geometries cannot represent:
#
# - WKT containing `EMPTY` (including nested, e.g.
#   `GEOMETRYCOLLECTION(POLYGON EMPTY, ...)`), because GI wrapper geometries
#   cannot be empty.
# - `GEOMETRYCOLLECTION`, because WellKnownGeometry mis-splits subgeometries
#   preceded by whitespace.
# - `LINEARRING`, which WellKnownGeometry does not know at all.
#
# The LibGEOS geometries are GeoInterface-compatible, so consumers (the
# RelateNG engine) access them via GI accessors like any other geometry.

import WellKnownGeometry
import GeoFormatTypes as GFT
import GeometryOps as GO
import LibGEOS as LG

function from_wkt(wkt::String)
    sanitized_wkt = join(strip.(split(wkt, "\n")), "")
    upper_wkt = uppercase(lstrip(sanitized_wkt))
    if occursin("EMPTY", upper_wkt) ||
            startswith(upper_wkt, "GEOMETRYCOLLECTION") ||
            startswith(upper_wkt, "LINEARRING")
        return LG.readgeom(sanitized_wkt)
    end
    geom = GFT.WellKnownText(GFT.Geom(), sanitized_wkt)
    return GO.tuples(geom)
end
