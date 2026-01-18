module GeometryOpsTestHelpersGeometryBasicsExt

using GeometryOpsTestHelpers
using GeoInterface
using GeometryBasics
import GeometryOps as GO

function __init__()
    # Register GeometryBasics in the test modules list
    push!(GeometryOpsTestHelpers.TEST_MODULES, GeometryBasics)
end

# Monkey-patch GeometryBasics to have correct methods.
# TODO: push this up to GB!

# TODO: remove when GB GI pr lands
function GeoInterface.convert(
    ::Type{GeometryBasics.LineString},
    ::GeoInterface.LinearRingTrait,
    geom
)
    return GeoInterface.convert(GeometryBasics.LineString, GeoInterface.LineStringTrait(), geom)
end
GeometryBasics.geointerface_geomtype(::GeoInterface.LinearRingTrait) = GeometryBasics.LineString

function GeoInterface.convert(::Type{GeometryBasics.Line}, type::GeoInterface.LineTrait, geom)
    g1, g2 = GeoInterface.getgeom(geom)
    x, y = GeoInterface.x(g1), GeoInterface.y(g1)
    if GeoInterface.is3d(geom)
        z = GeoInterface.z(g1)
        T = promote_type(typeof(x), typeof(y), typeof(z))
        return GeometryBasics.Line{3,T}(GeometryBasics.Point{3,T}(x, y, z), GeometryBasics.Point{3,T}(GeoInterface.x(g2), GeoInterface.y(g2), GeoInterface.z(g2)))
    else
        T = promote_type(typeof(x), typeof(y))
        return GeometryBasics.Line{2,T}(GeometryBasics.Point{2,T}(x, y), GeometryBasics.Point{2,T}(GeoInterface.x(g2), GeoInterface.y(g2)))
    end
end

# GeometryCollection interface - currently just a large Union
const _ALL_GB_GEOM_TYPES = Union{GeometryBasics.Point, GeometryBasics.LineString, GeometryBasics.Polygon, GeometryBasics.MultiPolygon, GeometryBasics.MultiLineString, GeometryBasics.MultiPoint}
GeometryBasics.geointerface_geomtype(::GeoInterface.GeometryCollectionTrait) = Vector{_ALL_GB_GEOM_TYPES}
function GeoInterface.convert(::Type{Vector{<: _ALL_GB_GEOM_TYPES}}, ::GeoInterface.GeometryCollectionTrait, geoms)
    return _ALL_GB_GEOM_TYPES[GeoInterface.convert(GeometryBasics, g) for g in GeoInterface.getgeom(geoms)]
end

function GeoInterface.convert(
    ::Type{GeometryBasics.LineString},
    type::GeoInterface.LineStringTrait,
    geom::GeoInterface.Wrappers.LinearRing{false, false, GO.StaticArrays.SVector{N, Tuple{Float64, Float64}}, Nothing, Nothing} where N
)
    return GeometryBasics.LineString(GeometryBasics.Point2{Float64}.(collect(geom.geom)))
end

end # module GeometryOpsTestHelpersGeometryBasicsExt
