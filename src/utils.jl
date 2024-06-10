# # Utility functions
_ismeasured(geom)::Bool = _ismeasured(GI.trait(geom), geom)
_ismeasured(::GI.AbstractGeometryTrait, geom)::Bool = GI.ismeasured(geom)
_ismeasured(::GI.FeatureTrait, feature)::Bool = _ismeasured(GI.geometry(feature))
_ismeasured(::GI.FeatureCollectionTrait, fc)::Bool = _ismeasured(GI.getfeature(fc, 1))
_ismeasured(::Nothing, geom)::Bool = _ismeasured(first(geom)) # Otherwise step into an itererable

_is3d(geom)::Bool = _is3d(GI.trait(geom), geom)
_is3d(::GI.AbstractGeometryTrait, geom)::Bool = GI.is3d(geom)
_is3d(::GI.FeatureTrait, feature)::Bool = _is3d(GI.geometry(feature))
_is3d(::GI.FeatureCollectionTrait, fc)::Bool = _is3d(GI.getfeature(fc, 1))
_is3d(::Nothing, geom)::Bool = _is3d(first(geom)) # Otherwise step into an itererable

_npoint(x) = _npoint(trait(x), x)
_npoint(::Nothing, xs::AbstractArray) = sum(_npoint, xs)
_npoint(::GI.FeatureCollectionTrait, fc) = sum(_npoint, GI.getfeature(fc))
_npoint(::GI.FeatureTrait, f) = _npoint(GI.geometry(f))
_npoint(::GI.AbstractGeometryTrait, x) = GI.npoint(trait(x), x)

_nedge(x) = _nedge(trait(x), x)
_nedge(::Nothing, xs::AbstractArray) = sum(_nedge, xs)
_nedge(::GI.FeatureCollectionTrait, fc) = sum(_nedge, GI.getfeature(fc))
_nedge(::GI.FeatureTrait, f) = _nedge(GI.geometry(f))
function _nedge(::GI.AbstractGeometryTrait, x)
    n = 0
    for g in GI.getgeom(x)
        n += _nedge(g)
    end
    return n
end
_nedge(::GI.AbstractCurveTrait, x) = GI.npoint(x) - 1
_nedge(::GI.PointTrait, x) = error("Cant get edges from points")


"""
    polygon_to_line(poly::Polygon)

Converts a Polygon to LineString or MultiLineString

# Examples

```jldoctest
import GeometryOps as GO, GeoInterface as GI

poly = GI.Polygon([[(-2.275543, 53.464547), (-2.275543, 53.489271), (-2.215118, 53.489271), (-2.215118, 53.464547), (-2.275543, 53.464547)]])
GO.polygon_to_line(poly)
# output
GeoInterface.Wrappers.LineString{false, false, Vector{Tuple{Float64, Float64}}, Nothing, Nothing}([(-2.275543, 53.464547), (-2.275543, 53.489271), (-2.215118, 53.489271), (-2.215118, 53.464547), (-2.275543, 53.464547)], nothing, nothing)
```
"""
function polygon_to_line(poly)
    @assert GI.trait(poly) isa PolygonTrait
    GI.ngeom(poly) > 1 && return GI.MultiLineString(collect(GI.getgeom(poly)))
    return GI.LineString(collect(GI.getgeom(GI.getgeom(poly, 1))))
end


"""
    to_edges()

Convert any geometry or collection of geometries into a flat 
vector of `Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}` edges.
"""
function to_edges(x, ::Type{T} = Float64) where T
    edges = Vector{TupleEdge{T}}(undef, _nedge(x))
    _to_edges!(edges, x, 1)
    return edges
end

_to_edges!(edges::Vector, x, n) = _to_edges!(edges, trait(x), x, n)
function _to_edges!(edges::Vector, ::GI.FeatureCollectionTrait, fc, n)
    for f in GI.getfeature(fc)
        n = _to_edges!(edges, f, n)
    end
end
_to_edges!(edges::Vector, ::GI.FeatureTrait, f, n) = _to_edges!(edges, GI.geometry(f), n)
function _to_edges!(edges::Vector, ::GI.AbstractGeometryTrait, fc, n)
    for f in GI.getgeom(fc)
        n = _to_edges!(edges, f, n)
    end
end
function _to_edges!(edges::Vector, ::GI.AbstractCurveTrait, geom, n)
    p1 = GI.getpoint(geom, 1) 
    p1x, p1y = GI.x(p1), GI.y(p1)
    for i in 2:GI.npoint(geom)
        p2 = GI.getpoint(geom, i)
        p2x, p2y = GI.x(p2), GI.y(p2)
        edges[n] = (p1x, p1y), (p2x, p2y)
        p1x, p1y = p2x, p2y
        n += 1
    end
    return n
end

_tuple_point(p) = GI.x(p), GI.y(p)
_tuple_point(p, ::Type{T}) where T = T(GI.x(p)), T(GI.y(p))

function to_extent(edges::Vector{<:Edge})
    x, y = extrema(first, edges)
    Extents.Extent(X=x, Y=y)
end

function to_points(x, ::Type{T} = Float64) where T
    points = Vector{TuplePoint{T}}(undef, _npoint(x))
    _to_points!(points, x, 1)
    return points
end

_to_points!(points::Vector, x, n) = _to_points!(points, trait(x), x, n)
function _to_points!(points::Vector, ::FeatureCollectionTrait, fc, n)
    for f in GI.getfeature(fc)
        n = _to_points!(points, f, n)
    end
end
_to_points!(points::Vector, ::FeatureTrait, f, n) = _to_points!(points, GI.geometry(f), n)
function _to_points!(points::Vector, ::AbstractGeometryTrait, fc, n)
    for f in GI.getgeom(fc)
        n = _to_points!(points, f, n)
    end
end
function _to_points!(points::Vector, ::Union{AbstractCurveTrait,MultiPointTrait}, geom, n)
    n = 0
    for p in GI.getpoint(geom)
        n += 1
        points[n] = _tuple_point(p)
    end
    return n
end

function _point_in_extent(p, extent::Extents.Extent)
    (x1, x2), (y1, y2) = extent.X, extent.Y
    return x1 ≤ GI.x(p) ≤ x2 && y1 ≤ GI.y(p) ≤ y2
end

_get_point_type(::Type{T}) where T = SVPoint_2D{T}

_sv_point(p, ::Type{T}) where T = SVPoint_2D{T}(_tuple_point(p))
# Get type of polygons that will be made
# TODO: Increase type options
_get_poly_type(::Type{T}) where T =
    GI.Polygon{false, false, Vector{GI.LinearRing{false, false, Vector{_get_point_type(T)}, Nothing, Nothing}}, Nothing, Nothing}

