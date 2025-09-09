# # GEOS Buffer

#=
## What is GEOS buffer?

The GEOS buffer function creates a new geometry that represents all points within
a specified distance of the input geometry. This is useful for creating zones of
influence, safety margins, or areas of interest around geometric features.

For example, creating a buffer around a point:

```@example buffer
import GeometryOps as GO
import GeoInterface as GI
using Makie
using CairoMakie

# Create a point and buffer it
point = GI.Point(0, 0)
buffered = GO.buffer(GO.GEOS(), point, 1.0)  # 1 unit buffer
```

## Implementation

The implementation uses GEOS's bufferWithStyle function through LibGEOS.jl.
It supports various buffer styles:
- Cap styles: round, flat, square
- Join styles: round, mitre, bevel
- Mitre limit for sharp corners

Key features:
- Configurable number of segments for curved parts (quadsegs)
- Support for different end cap styles
- Customizable join styles for corners
- Adjustable mitre limit for sharp angles
- Preserves CRS information
- Works with any GeoInterface-compatible geometry

The function handles the conversion between GeometryOps and GEOS geometries
automatically, and wraps the result back in a GeoInterface-compatible format.
=#

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

function GO.buffer(alg::GEOS, geometry, distance; calc_extent = true, kwargs...)
    # The reason we use apply here is so that this also works with featurecollections,
    # tables, vectors of geometries, etc!
    return apply(TraitTarget{GI.AbstractGeometryTrait}(), geometry; kwargs...) do geom
        newgeom = LG.bufferWithStyle(
            GI.convert(LG, geom), distance; 
            quadsegs = get(alg, :quadsegs, 8),
            endCapStyle = to_cap_style(get(alg, :endCapStyle, :round)),
            joinStyle = to_join_style(get(alg, :joinStyle, :round)),
            mitreLimit = get(alg, :mitreLimit, 5.0),
        )
        return _wrap(newgeom; crs = GI.crs(geom), calc_extent)
    end
end