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
GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}[GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([4.5, 3.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([6.5, 5.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([8.5, 7.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([4.5, 3.5], nothing)], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}[GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([6.5, 5.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([8.5, 7.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([9.5, 8.5], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([6.5, 5.5], nothing)], nothing, nothing)], nothing, nothing)
```

With Rotations.jl you need to actuall multiply the Rotation
by the `SVector` point, which is easy using an anonymous function.

```julia
julia> using Rotations

julia> GO.transform(p -> one(RotMatrix{2}) * p, geom)
GeoInterface.Wrappers.Polygon{false, false, Vector{GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}[GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}[GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([1.0, 2.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([3.0, 4.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([5.0, 6.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([1.0, 2.0], nothing)], nothing, nothing), GeoInterface.Wrappers.LinearRing{false, false, Vector{GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}}, Nothing, Nothing}(GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}[GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([3.0, 4.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([5.0, 6.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([6.0, 7.0], nothing), GeoInterface.Wrappers.Point{false, false, StaticArraysCore.SVector{2, Float64}, Nothing}([3.0, 4.0], nothing)], nothing, nothing)], nothing, nothing)
```
"""
function transform(f, geom, ::Type{T} = Float64; kw...) where T
    if _ismeasured(geom)
        return apply(PointTrait(), geom; kw...) do p
            GI.Point(T.(f(SA.SVector{4}(GI.x(p), GI.y(p), GI.z(p), GI.m(p)))))
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            GI.Point(T.(f(SA.SVector{3}((GI.x(p), GI.y(p), GI.z(p))))))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            GI.Point(T.(f(SA.SVector{2}((GI.x(p), GI.y(p))))))
        end
    end
end
