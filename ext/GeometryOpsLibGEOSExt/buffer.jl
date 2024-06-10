const _GEOS_CAPSTYLE_LOOKUP = Dict{Symbol, LG.GEOSBufCapStyles}(
    :round => LG.GEOSBUF_CAP_ROUND,
    :flat => LG.GEOSBUF_CAP_FLAT,
    :square => LG.GEOSBUF_CAP_SQUARE,
)

const _GEOS_JOINSTYLE_LOOKUP = Dict{Symbol, LG.GEOSBufJoinStyles}(
    :round => LG.GEOSBUF_JOIN_ROUND,
    :mitre => LG.GEOSBUF_JOIN_MITRE,
    :bevel => LG.GEOSBUF_JOIN_BEVEL,
)

to_cap_style(style::Symbol) = _GEOS_CAPSTYLE_LOOKUP[style]
to_cap_style(style::LG.GEOSBufCapStyles) = style
to_cap_style(num::Integer) = num

to_join_style(style::Symbol) = _GEOS_JOINSTYLE_LOOKUP[style]
to_join_style(style::LG.GEOSBufJoinStyles) = style
to_join_style(num::Integer) = num

function GO.buffer(alg::GEOS, geometry, distance)
    # The reason we use apply here is so that this also works with featurecollections,
    # tables, vectors of geometries, etc!
    return apply(TraitTarget{GI.AbstractGeometryTrait}(), geometry) do geom
        LG.bufferWithStyle(
            GI.convert(LG, geom), distance; 
            quadsegs = get(alg, :quadsegs, 8),
            endCapStyle = to_cap_style(get(alg, :endCapStyle, :round)),
            joinStyle = to_join_style(get(alg, :joinStyle, :round)),
            mitreLimit = get(alg, :mitreLimit, 5.0),
        )
    end
end