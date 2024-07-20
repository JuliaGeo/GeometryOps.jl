# # Pointwise transformation

"""
    transform(f, obj)

Apply a function `f` to all the points in `obj`.

Points will be passed to `f` as an `SVector` to allow
using CoordinateTransformations.jl and Rotations.jl 
without hassle.

`SVector` is also a valid GeoInterface.jl point, so will
work in all GeoInterface.jl methods.

## Example

```julia
julia> import GeoInterface as GI

julia> import GeometryOps as GO

julia> geom = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])]);

julia> f = CoordinateTransformations.Translation(3.5, 1.5)
Translation(3.5, 1.5)

julia> GO.transform(f, geom)
GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Linea
rRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}(StaticArraysCo
re.SVector{2, Float64}[[4.5, 3.5], [6.5, 5.5], [8.5, 7.5], [4.5, 3.5]], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Float64}}, Nothing, Nothing}(StaticA
rraysCore.SVector{2, Float64}[[6.5, 5.5], [8.5, 7.5], [9.5, 8.5], [6.5, 5.5]], nothing, nothing)], nothing, nothing)
```

With Rotations.jl you need to actually multiply the Rotation
by the `SVector` point, which is easy using an anonymous function.

```julia
julia> using Rotations

julia> GO.transform(p -> one(RotMatrix{2}) * p, geom)
GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.LinearR
ing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}(StaticArraysCore.SVe
ctor{2, Int64}[[2, 1], [4, 3], [6, 5], [2, 1]], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{StaticArraysCore.SVector{2, Int64}}, Nothing, Nothing}(StaticArraysCore.SVector{2, Int64
}[[4, 3], [6, 5], [7, 6], [4, 3]], nothing, nothing)], nothing, nothing)
```
"""
function transform(f, geom; kw...) 
    if _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            f(StaticArrays.SVector{3}((GI.x(p), GI.y(p), GI.z(p))))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            f(StaticArrays.SVector{2}((GI.x(p), GI.y(p))))
        end
    end
end
