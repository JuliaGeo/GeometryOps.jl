#=
# Perimeter

Computes the perimeter of a geometry.

Perimeter is not defined for points and multipoints.

For linestrings and linearrings, it is the sum of the lengths of the segments.
For polygons, it is the sum of the lengths of the exterior and holes of the polygon.
For multipolygons, it is the sum of the perimeters of the polygons in the multipolygon.
For geometry collections, it is the sum of the perimeters of the geometries in the collection. 
The geometry collections cannot have point or multipoint geometries.

TODOs:
- Add support for point geometries
- Add support for a "true 3D" perimeter

```@example perimeter
import GeoInterface as GI, GeometryOps as GO

outer = GI.LinearRing([(0,0),(10,0),(10,10),(0,10),(0,0)])
hole1 = GI.LinearRing([(1,1),(1,2),(2,2),(2,1),(1,1)])
hole2 = GI.LinearRing([(5,5),(5,6),(6,6),(6,5),(5,5)])

p = GI.Polygon([outer, hole1, hole2])
mp = GI.MultiPolygon([
    p, 
    GO.transform(x -> x .+ 12, GI.Polygon([outer, hole1]))
])

(p, mp)
```
```@example perimeter
GO.perimeter(p) # should be 48
```
```@example perimeter
GO.perimeter(mp) # should be 92
```
=#

const _PERIMETER_TARGETS = TraitTarget{GI.AbstractCurveTrait}()

"""
    abstract type DistanceAlgorithm

An abstract supertype for a distance algorithm that computes the distance between two points.

Currently used in [`GO.perimeter`](@ref GeometryOps.perimeter), but should be used in GO.distance and other functions as well...

See also: [`LinearDistance`](@ref), [`GeodesicDistance`](@ref), [`RhumbDistance`](@ref).
"""
abstract type DistanceAlgorithm end

"""
    LinearDistance()

A linear distance algorithm that uses simple 2D Euclidean distance between points.
"""
struct LinearDistance <: DistanceAlgorithm end

"""
    GeodesicDistance()

A geodesic distance algorithm that uses the geodesic distance between points.

Requires the Proj.jl package to be loaded, uses Proj's GeographicLib.
"""
struct GeodesicDistance{T} <: DistanceAlgorithm 
    geodesic::T# ::Proj.geod_geodesic
end

"""
    RhumbDistance()

A rhumb distance algorithm that uses the rhumb distance between points.
"""
struct RhumbDistance <: DistanceAlgorithm end

"""
    perimeter([alg::DistanceAlgorithm], geom, [T::Type{<: Number} = Float64]; threaded = false)

Computes the perimeter of a geometry using the specified distance algorithm.  

Allowed distance algorithms are: [`LinearDistance`](@ref) (default), [`GeodesicDistance`](@ref), [`RhumbDistance`](@ref).  The latter two are not yet implemented.
"""
perimeter(geom, T::Type{_T} = Float64; threaded::Union{Bool, BoolsAsTypes} = _False()) where _T <: Number = perimeter(LinearDistance(), geom, _T; threaded = threaded)

function perimeter(alg::DistanceAlgorithm, geom, T::Type{_T} = Float64; threaded::Union{Bool, BoolsAsTypes} = _False()) where _T <: Number
    _threaded = _booltype(threaded)
    find_perimeter(geom) = _perimeter(alg, geom, T)
    return applyreduce(find_perimeter, +, _PERIMETER_TARGETS, geom; threaded = _threaded, init = zero(T))
end



_perimeter(alg::DistanceAlgorithm, geom, ::Type{T}) where T <: Number = _perimeter(GI.trait(geom), alg, geom, T)

function _perimeter(::GI.AbstractCurveTrait, alg::DistanceAlgorithm, geom, ::Type{T}) where T <: Number
    ret = T(0)
    prev = GI.getpoint(geom, 1)
    for point in GI.getpoint(geom)
        ret += point_distance(alg, prev, point, T)
        prev = point
    end
    return ret
end

point_distance(::LinearDistance, p1, p2, ::Type{T}) where T <: Number = T(hypot(GI.x(p2) - GI.x(p1), GI.y(p2) - GI.y(p1)))
point_distance(alg::DistanceAlgorithm, p1, p2, ::Type{T}) where T <: Number = error("Not implemented yet for alg $alg")

# point_distance(::RhumbDistance, p1, p2, ::Type{T}) where T <: Number = ...

# Add an error hint for Geodesicdistance if Proj is not loaded!
function _geodesic_distance_error_hinter(io, exc, argtypes, kwargs)
    if isnothing(Base.get_extension(GeometryOps, :GeometryOpsProjExt)) && exc.f == GeodesicDistance
        print(io, "\n\nThe `GeodesicDistance` method requires the Proj.jl package to be explicitly loaded.\n")
        print(io, "You can do this by simply typing ")
        printstyled(io, "using Proj"; color = :cyan, bold = true)
        println(io, " in your REPL, \nor otherwise loading Proj.jl via using or import.")
    end
end