# # Utility functions

_is3d(geom; kw...)::Bool = _is3d(GI.trait(geom), geom; kw...)
_is3d(::GI.AbstractGeometryTrait, geom; geometrycolumn = nothing)::Bool = GI.is3d(geom)
_is3d(::GI.FeatureTrait, feature; geometrycolumn = nothing)::Bool = _is3d(GI.geometry(feature))
_is3d(::GI.FeatureCollectionTrait, fc; geometrycolumn = nothing)::Bool = _is3d(GI.getfeature(fc, 1))
function _is3d(::Nothing, geom; geometrycolumn = nothing)::Bool
    if Tables.istable(geom)
        geometrycolumn = isnothing(geometrycolumn) ? GI.geometrycolumns(geom) : geometrycolumn isa Symbol ? (geometrycolumn,) : geometrycolumn
        # take the first geometry column
        # TODO: this is a bad guess - this should really be on the vector level somehow.  
        # Maybe a configurable applicator again....
        first_geom = if Tables.rowaccess(geom)
            Tables.getcolumn(first(Tables.rows(geom)), first(geometrycolumn))
        else # column access assumed
            first(Tables.getcolumn(geom, first(geometrycolumn)))
        end
        return _is3d(first_geom)
    else # assume iterable
        first_geom = first(geom)
        if GI.trait(first_geom) isa GI.AbstractTrait
            return _is3d(first_geom)
        else
            return false # couldn't figure it out!
        end
    end

end

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
# implementation for Polygon, MultiPolygon, MultiLineString, GeometryCollection
function eachedge(trait::GI.AbstractGeometryTrait, geom, ::Type{T}) where T
    return Iterators.flatten((eachedge(r, T) for r in flatten(GI.AbstractCurveTrait, geom)))
end
function eachedge(trait::GI.PointTrait, geom, ::Type{T}) where T
    throw(ArgumentError("Can't get edges from points, $geom was a PointTrait."))
end
function eachedge(trait::GI.MultiPointTrait, geom, ::Type{T}) where T
    throw(ArgumentError("Can't get edges from MultiPoint, $geom was a MultiPointTrait."))
end

"""
    to_edgelist(geom, [::Type{T}])

Convert a geometry into a vector of `GI.Line` objects with attached extents.
"""
to_edgelist(geom, ::Type{T} = Float64) where T = 
    [_lineedge(ps, T) for ps in eachedge(geom, T)]

"""
    to_edgelist(ext::E, geom, [::Type{T}])::(::Vector{GI.Line}, ::Vector{Int})

Filter the edges of `geom` for those that intersect 
`ext`, and return:
- a vector of `GI.Line` objects with attached extents, 
- a vector of indices into the original geometry.
"""
function to_edgelist(ext::E, geom, ::Type{T} = Float64) where {E<:Extents.Extent,T}
    edges_in = eachedge(geom, T)
    l1 = _lineedge(first(edges_in), T)
    edges_out = typeof(l1)[]
    indices = Int[]
    for (i, ps) in enumerate(edges_in) 
        l = _lineedge(ps, T)
        if Extents.intersects(ext, GI.extent(l))
            push!(edges_out, l)
            push!(indices, i)
        end
    end 
    return edges_out, indices
end

function _lineedge(ps::Tuple, ::Type{T}) where T
    l = GI.Line(StaticArrays.SVector{2,NTuple{2, T}}(ps))  # TODO: make this flexible in dimension
    e = GI.extent(l)
    return GI.Line(l.geom; extent=e)
end

"""
    lazy_edgelist(geom, [::Type{T}])

Return an iterator over `GI.Line` objects with attached extents.
"""
function lazy_edgelist(geom, ::Type{T} = Float64) where T
    (_lineedge(ps, T) for ps in eachedge(geom, T))
end

"""
    edge_extents(geom, [::Type{T}])

Return a vector of the extents of the edges (line segments) of `geom`.
"""
function edge_extents(geom, ::Type{T} = Float64) where T
    return [begin
        Extents.Extent(X=extrema(GI.x, edge), Y=extrema(GI.y, edge)) 
    end
    for edge in eachedge(geom, T)]
end

"""
    lazy_edge_extents(geom)

Return an iterator over the extents of the edges (line segments) of `geom`.  
This is lazy but nonallocating.
"""
function lazy_edge_extents(geom)
    return (begin
        Extents.Extent(X=extrema(GI.x, edge), Y=extrema(GI.y, edge)) 
    end
    for edge in eachedge(geom, Float64))
end


# Extent to polygon

"""
    extent_to_polygon(ext::Extents.Extent)

Convert an extent to a polygon.

# Examples

```jldoctest
import GeometryOps as GO, Extents
    
ext = Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0))
GO.extent_to_polygon(ext)
# output
GeoInterface.Wrappers.Polygon{false, false}([GeoInterface.Wrappers.LinearRing([(1.0, 1.0), … (3) … , (1.0, 1.0)])])
```
"""
function extent_to_polygon(ext::Extents.Extent{(:X, :Y)})
    x1, x2 = ext.X
    y1, y2 = ext.Y
    return GI.Polygon(StaticArrays.@SVector[GI.LinearRing(StaticArrays.@SVector[(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])
end


function extent_to_polygon(ext::Extents.Extent{(:Y, :X)})
    x1, x2 = ext.X
    y1, y2 = ext.Y
    return GI.Polygon(StaticArrays.@SVector[GI.LinearRing(StaticArrays.@SVector[(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])
end

# This will accept table rows etc
# TODO can rows have metadata to detectg the geometry column name?
function _geometry_or_error(g2)
    hasproperty(g2, :geometry) || throw(ArgumentError("Objects that return no geometry or feature traits must at least have a :geometry property"))
    return g2.geometry
end


