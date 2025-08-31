# # Arclength

export arclength_to_point, point_at_arclength

#=
## What is arclength functionality?

Arclength functionality provides two key operations:
1. `arclength_to_point(manifold, linestring, point)` - calculates the cumulative 
   distance along a linestring from its start to a specified point on the line
2. `point_at_arclength(manifold, linestring, distance)` - finds the point at a 
   specified distance along the linestring from its start

These functions are useful for:
- Parameterizing curves by arc length
- Finding positions along routes or paths
- Interpolating along geometric curves
- Measuring progress along linear features

Both functions support multiple manifolds:
- `Planar()` - uses Euclidean distance calculations
- `Geodesic()` - uses geodesic distance calculations for accurate Earth-surface measurements

## Examples

```@example arclength
import GeometryOps as GO, GeoInterface as GI

# Create a simple linestring
line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0)])

# Find the distance to a point on the line
distance_to_point = GO.arclength_to_point(GO.Planar(), line, (1.0, 0.5))

# Find a point at a specific distance along the line
point_at_distance = GO.point_at_arclength(GO.Planar(), line, 1.5)
```

For geographic coordinates, use the Geodesic manifold:

```@example arclength
using Proj # required for Geodesic calculations
geo_line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
distance_geo = GO.arclength_to_point(GO.Geodesic(), geo_line, (0.5, 0.0))
point_geo = GO.point_at_arclength(GO.Geodesic(), geo_line, 50000) # 50km
```

## Implementation
=#

"""
    arclength_to_point([method = Planar()], linestring, point; threaded = false)

Calculate the cumulative distance along a linestring from its start to a 
specified point. The point should lie on the linestring.

## Arguments
- `method::Manifold = Planar()`: The manifold to use for distance calculations.
  - `Planar()` uses Euclidean distance
  - `Geodesic()` uses geodesic distance calculations
- `linestring`: A LineString or LinearRing geometry
- `point`: The target point on the linestring

Returns the cumulative distance from the start of the linestring to the point.
If the point is not on the linestring, returns the distance to the closest point
on the linestring.

## Example
```julia
import GeometryOps as GO, GeoInterface as GI

line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
distance = GO.arclength_to_point(line, (1.0, 0.5))
```
"""
function arclength_to_point(linestring, point; threaded::Union{Bool, BoolsAsTypes} = False())
    return arclength_to_point(Planar(), linestring, point; threaded = booltype(threaded))
end

function arclength_to_point(method::Manifold, linestring, point; threaded::Union{Bool, BoolsAsTypes} = False())
    return _arclength_to_point(method, linestring, point, GI.trait(linestring))
end

"""
    point_at_arclength([method = Planar()], linestring, distance; threaded = false)

Find the point at a specified distance along a linestring from its start.

## Arguments
- `method::Manifold = Planar()`: The manifold to use for distance calculations.
  - `Planar()` uses Euclidean distance
  - `Geodesic()` uses geodesic distance calculations  
- `linestring`: A LineString or LinearRing geometry
- `distance`: The target distance along the linestring

Returns the point at the specified distance. If the distance exceeds the total
length of the linestring, returns the endpoint.

## Example
```julia
import GeometryOps as GO, GeoInterface as GI

line = GI.LineString([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
point = GO.point_at_arclength(line, 1.5)
```
"""
function point_at_arclength(linestring, distance; threaded::Union{Bool, BoolsAsTypes} = False())
    return point_at_arclength(Planar(), linestring, distance; threaded = booltype(threaded))
end

function point_at_arclength(method::Manifold, linestring, distance; threaded::Union{Bool, BoolsAsTypes} = False())
    return _point_at_arclength(method, linestring, distance, GI.trait(linestring))
end

# Implementation for LineString and LinearRing
function _arclength_to_point(method::Union{Planar, Spherical}, linestring, target_point, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    cumulative_distance = 0.0
    closest_distance = Inf
    result_distance = 0.0
    
    if GI.npoint(linestring) < 2
        return 0.0
    end
    
    prev_point = GI.getpoint(linestring, 1)
    
    for i in 2:GI.npoint(linestring)
        curr_point = GI.getpoint(linestring, i)
        
        # Calculate distance between consecutive points
        segment_length = _point_distance(method, prev_point, curr_point)
        
        # Find closest point on this segment to target
        closest_point_on_segment, t = _closest_point_on_segment(method, prev_point, curr_point, target_point)
        distance_to_segment = _point_distance(method, target_point, closest_point_on_segment)
        
        # If this is the closest segment so far
        if distance_to_segment < closest_distance
            closest_distance = distance_to_segment
            # Calculate distance to the closest point on this segment
            result_distance = cumulative_distance + t * segment_length
        end
        
        cumulative_distance += segment_length
        prev_point = curr_point
    end
    
    return result_distance
end

# Geodesic implementation is in the Proj extension
function _arclength_to_point(method::Geodesic, linestring, target_point, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    error("Geodesic arclength calculation requires Proj.jl to be loaded. Please run `using Proj`.")
end

function _point_at_arclength(method::Union{Planar, Spherical}, linestring, target_distance, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    if GI.npoint(linestring) < 2
        return GI.npoint(linestring) > 0 ? GI.getpoint(linestring, 1) : nothing
    end
    
    if target_distance <= 0
        return GI.getpoint(linestring, 1)
    end
    
    cumulative_distance = 0.0
    prev_point = GI.getpoint(linestring, 1)
    
    for i in 2:GI.npoint(linestring)
        curr_point = GI.getpoint(linestring, i)
        segment_length = _point_distance(method, prev_point, curr_point)
        
        if cumulative_distance + segment_length >= target_distance
            # Target distance is within this segment
            remaining_distance = target_distance - cumulative_distance
            t = remaining_distance / segment_length
            return _interpolate_point(method, prev_point, curr_point, t)
        end
        
        cumulative_distance += segment_length
        prev_point = curr_point
    end
    
    # Distance exceeds total length, return last point
    return GI.getpoint(linestring, GI.npoint(linestring))
end

# Geodesic implementation is in the Proj extension
function _point_at_arclength(method::Geodesic, linestring, target_distance, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    error("Geodesic point at arclength calculation requires Proj.jl to be loaded. Please run `using Proj`.")
end

# Helper functions for distance calculations
function _point_distance(::Planar, p1, p2)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    return hypot(x2 - x1, y2 - y1)
end

# Geodesic version is defined in the Proj extension
function _point_distance(::Geodesic, p1, p2, proj_geodesic)
    error("Geodesic distance calculation requires Proj.jl to be loaded. Please run `using Proj`.")
end

# Find the closest point on a line segment to a target point
function _closest_point_on_segment(::Planar, p1, p2, target)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    tx, ty = GI.x(target), GI.y(target)
    
    # Vector from p1 to p2
    dx = x2 - x1
    dy = y2 - y1
    
    # Vector from p1 to target
    px = tx - x1
    py = ty - y1
    
    # Project target onto line segment
    segment_length_sq = dx * dx + dy * dy
    
    if segment_length_sq == 0
        # Degenerate segment
        return p1, 0.0
    end
    
    t = (px * dx + py * dy) / segment_length_sq
    t = clamp(t, 0.0, 1.0)  # Clamp to segment bounds
    
    closest_x = x1 + t * dx
    closest_y = y1 + t * dy
    
    return (closest_x, closest_y), t
end

# Geodesic version is defined in the Proj extension
function _closest_point_on_segment(::Geodesic, p1, p2, target, proj_geodesic)
    error("Geodesic closest point calculation requires Proj.jl to be loaded. Please run `using Proj`.")
end

# Interpolate between two points
function _interpolate_point(::Planar, p1, p2, t)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    
    return (x, y)
end

# Geodesic version is defined in the Proj extension
function _interpolate_point(::Geodesic, p1, p2, t, proj_geodesic)
    error("Geodesic interpolation requires Proj.jl to be loaded. Please run `using Proj`.")
end