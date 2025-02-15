# # Utility functions

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
GeoInterface.Wrappers.LineString{false, false}([(-2.275543, 53.464547), … (3) … , (-2.275543, 53.464547)])
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
    edges = Vector{Edge{T}}(undef, _nedge(x))
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

function to_extent(edges::Vector{Edge})
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

#=
# `eachedge`, `to_edgelist`

These functions are used to decompose geometries into lists of edges.
Currently they only work on linear rings.
=#

"""
    eachedge(geom, [::Type{T}])

Decompose a geometry into a list of edges.
Currently only works for LineString and LinearRing.

Returns some iterator, which yields tuples of points.  Each tuple is an edge.

It goes `(p1, p2), (p2, p3), (p3, p4), ...` etc.
"""
eachedge(geom) = eachedge(GI.trait(geom), geom, Float64)

function eachedge(geom, ::Type{T}) where T
    eachedge(GI.trait(geom), geom, T)
end

# implementation for LineString and LinearRing
function eachedge(trait::GI.AbstractCurveTrait, geom, ::Type{T}) where T
    return (_tuple_point.((GI.getpoint(geom, i), GI.getpoint(geom, i+1)), T) for i in 1:GI.npoint(geom)-1)
end

"""
    to_edgelist(geom, [::Type{T}])

Convert a geometry into a vector of `GI.Line` objects with attached extents.
"""
to_edgelist(geom, ::Type{T}) where T = 
    [_lineedge(ps, T) for ps in eachedge(geom, T)]


function to_edgelist(ext::E, geom, ::Type{T}) where {E<:Extent,T}
    indices = Int[]
    edges_in = eachedge(geom, T)
    l1 = _lineedge(first(edges_in), T)
    edges_out = typeof(l1)[]
    for (i, ps) in enumerate(edges_in) 
        l = _lineedge(ps, T)
        if _intersects(ext, l, T)
            push!(edges_out, l)
            push!(indices, i)
        end
    end 
    return edges_out, indices
end

function _lineedge(ps::Tuple, ::Type{T}) where T
    l = GI.Line(SVector{2,NTuple{2, T}}(ps))  # TODO: make this flexible in dimension
    e = GI.extent(l)
    return GI.Line(l.geom; extent=e)
end

function _intersects(ext::Extent, l::GI.Line, ::Type{T}) where T
    p1, p2 = GI.getpoint(l)
    # Check if the points are in the extent
    ext.X[1] <= p1[1] <= ext.X[2] && ext.Y[1] <= p1[2] <= ext.Y[2] || 
    ext.X[1] <= p2[1] <= ext.X[2] && ext.Y[1] <= p2[2] <= ext.Y[2] && return true
    Extents.intersects(ext, GI.extent(l)) && return true
    # Otherwise check if the line intersects the extent square
    a = GI.Line(SVector{2,NTuple{2,T}}((ext.X[1], ext.Y[1]), (ext.X[1], ext.Y[2])))
    b = GI.Line(SVector{2,NTuple{2,T}}((ext.X[2], ext.Y[1]), (ext.X[2], ext.Y[2])))
    c = GI.Line(SVector{2,NTuple{2,T}}((ext.X[1], ext.Y[1]), (ext.X[2], ext.Y[1])))
    d = GI.Line(SVector{2,NTuple{2,T}}((ext.X[1], ext.Y[2]), (ext.X[2], ext.Y[2])))
    return intersects(a, l) || intersects(b, l) || intersects(c, l) || intersects(d, l)
end