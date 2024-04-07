function GO.buffer(::GEOS, geom, distance)
    return LG.buffer(GI.convert(LG, geom), distance)
end

function GO.buffer(::GEOS{(:cap_style, :join_style)}, geom, distance)
    return LG.bufferWithStyle(
        GI.convert(LG, geom), distance; 
        endCapStyle = alg.params.cap_style, 
        joinStyle = alg.params.join_style
    )
end