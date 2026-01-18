# Rotating geometry

This tutorial demonstrates how to rotate geometries, particularly polygons, using GeometryOps.jl.

There are several approaches to rotating geometry, depending on your needs:

1. **Simple rotation using the `rotate` convenience function**
2. Simple rotation around the origin
3. Rotation around a polygon's centroid  
4. Using external transformation libraries

All rotation approaches ultimately use the [`transform`](@ref) function which applies a transformation to every point in a geometry.

## Using the rotate convenience function

The simplest way to rotate a geometry is using the `rotate` function:

```@example rotation
import GeometryOps as GO
import GeoInterface as GI
using CairoMakie

# Create a simple square polygon  
square = GI.Polygon([[(1, 1), (2, 1), (2, 2), (1, 2), (1, 1)]])

# Rotate by 45 degrees (π/4 radians) around the centroid (default)
rotated_square = GO.rotate(square, π/4)

# Rotate around a specific point (e.g., origin)
rotated_around_origin = GO.rotate(square, π/4; origin = (0, 0))

# Visualize the results
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Rotation using rotate() function")
poly!(ax, square, color = :blue, alpha = 0.3, label = "Original")
poly!(ax, rotated_square, color = :red, alpha = 0.5, label = "Rotated around centroid")
poly!(ax, rotated_around_origin, color = :green, alpha = 0.5, label = "Rotated around origin")
axislegend(ax)
fig
```

## Simple rotation around origin

The simplest approach is to rotate around the origin (0, 0) using a 2D rotation matrix.

```@example rotation
import GeometryOps as GO
import GeoInterface as GI
using CoordinateTransformations, StaticArrays
using CairoMakie

# Create a simple square polygon
square = GI.Polygon([[(1, 1), (2, 1), (2, 2), (1, 2), (1, 1)]])

# Create a 2D rotation matrix function (rotate by angle in radians)
rotate2d(angle) = StaticArrays.SMatrix{2,2}(cos(angle), sin(angle), -sin(angle), cos(angle))

# Rotate by 45 degrees (π/4 radians)
rotated_square = GO.transform(p -> rotate2d(π/4) * p, square)

# Visualize the original and rotated polygons
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Rotation around origin")
poly!(ax, [square, rotated_square], color = [:blue, :red], alpha = 0.5)
fig
```

## Rotation around polygon centroid

Often you want to rotate a polygon around its own center rather than the origin. This requires:

1. Translating the polygon so its centroid is at the origin
2. Applying the rotation
3. Translating back to the original centroid position

```@example rotation  
# Calculate the centroid of our square
centroid_point = GO.centroid(square)

# Create a rotation transformation around the centroid
function rotate_around_centroid(geom, angle)
    center = GO.centroid(geom)
    rotation_matrix = rotate2d(angle)
    
    return GO.transform(geom) do point
        # Translate to origin, rotate, then translate back
        rotated_point = rotation_matrix * (point .- center)
        return rotated_point .+ center
    end
end

# Rotate the square 90 degrees around its centroid
rotated_around_center = rotate_around_centroid(square, π/2)

# Visualize
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Rotation around centroid")
poly!(ax, [square, rotated_around_center], color = [:blue, :red], alpha = 0.5)
# Mark the centroid
scatter!(ax, [centroid_point[1]], [centroid_point[2]], color = :black, marker = :x, markersize = 15)
fig
```

## Using CoordinateTransformations.jl

For more complex transformations, you can use the CoordinateTransformations.jl library, which provides composable transformations.

```@example rotation
using CoordinateTransformations

# Create a more complex polygon
polygon = GI.Polygon([[(0, 0), (3, 0), (3, 2), (1, 2), (1, 1), (0, 1), (0, 0)]])

# Get the centroid
center = GO.centroid(polygon)

# Create a composite transformation: translate to origin, rotate, translate back  
rotation_transform = Translation(center) ∘ LinearMap(rotate2d(π/3)) ∘ Translation(-center[1], -center[2])

# Apply the transformation
rotated_polygon = GO.transform(rotation_transform, polygon)

# Visualize multiple rotations
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Multiple rotations using CoordinateTransformations")

# Show original
poly!(ax, polygon, color = :blue, alpha = 0.3, label = "Original")

# Show rotations at different angles
angles = [π/6, π/3, π/2, 2π/3]
colors = [:red, :green, :orange, :purple]

for (angle, color) in zip(angles, colors)
    transform = Translation(center) ∘ LinearMap(rotate2d(angle)) ∘ Translation(-center[1], -center[2])
    rotated = GO.transform(transform, polygon)
    poly!(ax, rotated, color = color, alpha = 0.5, label = "$(round(rad2deg(angle)))°")
end

axislegend(ax)
fig
```

## Using Rotations.jl

For 3D rotations or more sophisticated rotation operations, you can use Rotations.jl. However, since Rotations.jl objects are not directly callable, you need to wrap them:

```@example rotation
using Rotations

# For 2D rotation using Rotations.jl, we need to work with the rotation matrix
# and extract the 2D part
polygon_2d = GI.Polygon([[(1, 0), (2, 0), (2, 1), (1, 1), (1, 0)]])

# Create a Z-axis rotation (for 2D rotation in the XY plane)
function rotate_with_rotations_jl(geom, angle_degrees)
    # Create 3D rotation around Z axis
    rotation_3d = RotZ(deg2rad(angle_degrees))
    
    # Extract 2D rotation matrix (top-left 2x2 submatrix)
    rotation_2d = rotation_3d[1:2, 1:2]
    
    # Apply to geometry
    return GO.transform(p -> rotation_2d * p, geom)
end

# Rotate by different angles
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Rotation using Rotations.jl")

poly!(ax, polygon_2d, color = :blue, alpha = 0.3, label = "Original")

for (angle, color) in zip([45, 90, 135], [:red, :green, :orange])
    rotated = rotate_with_rotations_jl(polygon_2d, angle)
    poly!(ax, rotated, color = color, alpha = 0.5, label = "$(angle)°")
end

axislegend(ax)
fig
```

## Performance considerations

When rotating many geometries or performing many rotations, consider:

1. Pre-computing the rotation matrix
2. Using StaticArrays for better performance
3. Using the `threaded=true` option in transform for large geometries

```@example rotation
# Example of efficient batch rotation
polygons = [GI.Polygon([[(i, 0), (i+1, 0), (i+1, 1), (i, 1), (i, 0)]]) for i in 1:5]

# Pre-compute rotation matrix
rotation_mat = rotate2d(π/4)

# Rotate all polygons efficiently
rotated_polygons = [GO.transform(p -> rotation_mat * p, poly) for poly in polygons]

# For very large geometries, use threading
# large_polygon = ... # some large polygon
# rotated_large = GO.transform(p -> rotation_mat * p, large_polygon; threaded=true)

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect(), title = "Batch rotation")
poly!(ax, polygons, color = :blue, alpha = 0.3)
poly!(ax, rotated_polygons, color = :red, alpha = 0.5)
fig
```

## Summary

GeometryOps.jl provides flexible geometry rotation through:

- **`rotate(geom, angle)`** - Simple convenience function for common rotation tasks
- **`rotate(geom, angle; origin=point)`** - Rotate around a specific point  
- **`transform`** with rotation matrices for advanced control
- **CoordinateTransformations.jl** integration for complex, composable transformations
- **Rotations.jl** integration for advanced 3D rotations

The key insight is that both `rotate` and `transform` work with any function that takes a point (as an SVector) and returns a transformed point, making the system highly flexible for various geometric transformations.

For most use cases, `GO.rotate(geom, angle)` provides the simplest interface, while `transform` offers maximum flexibility for custom transformations.