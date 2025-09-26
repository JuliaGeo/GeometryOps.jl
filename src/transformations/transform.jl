# # Transform
#=

Transform a geometry by applying some function to all the points,
which are provided as StaticVectors (which math like `+`, `*`, etc. works on).

You can pass any function or object that takes an `SVector` and returns some GeoInterface-compatible point.

## Example

```@example transformations
import GeoInterface as GI, GeometryOps as GO

geom = GI.Polygon([[(0,0), (1,0), (1,1), (0,1), (0,0)]])
geom2 = GO.transform(p -> p .+ (1, 2), geom)

using CairoMakie, GeoInterfaceMakie
poly([geom, geom2]; color = [:steelblue, :orange])
``` 

This uses [`apply`](@ref), so will work with any geometry, vector of geometries, table, etc.
=#

export transform, rotate

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

## Rotation Examples

For rotating geometry, you can use various approaches:

### Simple rotation using the rotate convenience function
```julia
julia> # Rotate by 45 degrees around the centroid (default)
julia> GO.rotate(geom, π/4)

julia> # Rotate around a specific point
julia> GO.rotate(geom, π/4; origin = (0, 0))
```

### Simple 2D rotation around origin
```julia
julia> using StaticArrays

julia> # Rotate by 45 degrees (π/4 radians)
julia> rotation_matrix = @SMatrix [cos(π/4) -sin(π/4); sin(π/4) cos(π/4)]

julia> GO.transform(p -> rotation_matrix * p, geom)
```

### Rotation around geometry centroid
```julia
julia> # Rotate around the polygon's centroid
julia> center = GO.centroid(geom)

julia> GO.transform(geom) do p
           # Translate to origin, rotate, then translate back
           rotated = rotation_matrix * (p .- center)  
           return rotated .+ center
       end
```

### Using CoordinateTransformations.jl for complex rotations
```julia
julia> using CoordinateTransformations

julia> center = GO.centroid(geom)

julia> # Compose transformations: translate, rotate, translate back
julia> rotation_transform = Translation(center) ∘ LinearMap(rotation_matrix) ∘ Translation(-center[1], -center[2])

julia> GO.transform(rotation_transform, geom)
```

### Using Rotations.jl for 2D rotation
```julia
julia> using Rotations

julia> # For 2D rotation, extract 2x2 submatrix from 3D rotation
julia> rotation_2d = RotZ(π/4)[1:2, 1:2]

julia> GO.transform(p -> rotation_2d * p, geom)
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

"""
    rotate(geom, angle::Real; origin = nothing)

Rotate a geometry by `angle` (in radians) around a point.

If `origin` is not provided, the geometry is rotated around its centroid.
If `origin` is provided as a point-like object (e.g., tuple or Point), 
the geometry is rotated around that point.

## Examples

```julia
# Rotate a square by 45 degrees around its centroid
square = GI.Polygon([[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]])
rotated = GO.rotate(square, π/4)

# Rotate around a specific point
rotated_around_origin = GO.rotate(square, π/4; origin = (0, 0))

# Rotate by 90 degrees (π/2 radians)
rotated_90 = GO.rotate(square, π/2)
```
"""
function rotate(geom, angle::Real; origin = nothing)
    # Create 2D rotation matrix
    cos_a, sin_a = cos(angle), sin(angle)
    rotation_matrix = StaticArrays.@SMatrix [cos_a -sin_a; sin_a cos_a]
    
    if origin === nothing
        # Rotate around centroid
        center = centroid(geom)
        return transform(geom) do p
            rotated = rotation_matrix * (p .- center)
            return rotated .+ center
        end
    else
        # Rotate around specified origin point
        origin_point = StaticArrays.SVector{2}(origin)
        return transform(geom) do p
            rotated = rotation_matrix * (p .- origin_point)
            return rotated .+ origin_point
        end
    end
end
