# # Utility functions

_is3d(geom) = _is3d(GI.trait(geom), geom)
_is3d(::GI.AbstractGeometryTrait, geom) = GI.is3d(geom)
_is3d(::GI.FeatureTrait, feature) = _is3d(GI.geometry(feature))
_is3d(::GI.FeatureCollectionTrait, fc) = _is3d(GI.getfeature(fc, 1))
_is3d(::Nothing, geom) = _is3d(first(geom)) # Otherwise step into an itererable


"""
    polygon_to_line(poly::Polygon)

Converts a Polygon to LineString or MultiLineString

# Examples
```jldoctest
julia> poly = Polygon([[[-2.275543, 53.464547],[-2.275543, 53.489271],[-2.215118, 53.489271],[-2.215118, 53.464547],[-2.275543, 53.464547]]])
Polygon(Array{Array{Float64,1},1}[[[-2.27554, 53.4645], [-2.27554, 53.4893], [-2.21512, 53.4893], [-2.21512, 53.4645], [-2.27554, 53.4645]]])

julia> polygon_to_line(poly)
LineString(Array{Float64,1}[[-2.27554, 53.4645], [-2.27554, 53.4893], [-2.21512, 53.4893], [-2.21512, 53.4645], [-2.27554, 53.4645]])
```
"""
function polygon_to_line(poly)
    @assert GI.trait(poly) isa PolygonTrait
    GI.ngeom(poly) > 1 && return GI.MultiLineString(collect(GI.getgeom(poly)))
    return GI.LineString(collect(GI.getgeom(GI.getgeom(poly, 1))))
end
