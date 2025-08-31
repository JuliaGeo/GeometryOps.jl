#=

Geodesic arclength functionality via PROJ.

```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
GeometryOps.arclength_to_point
GeometryOps.point_at_arclength
```

Implementation

The implementation uses PROJ's geodesic calculations to:
1. Compute the geodesic distance between points accurately on the ellipsoid
2. Find closest points on geodesic segments to target points
3. Interpolate positions along geodesic lines at specified distances

Key features:
- Uses PROJ's geod_geodesic for accurate ellipsoidal calculations
- Configurable equatorial radius and flattening via Geodesic manifold
- Thread-safe implementation  
- Supports both LineString and LinearRing geometries

The function creates geodesic lines between each pair of points and
calculates distances and interpolated positions along those geodesic paths.
=#

# This holds the `arclength` geodesic functionality.

import GeometryOps: _arclength_to_point, _point_at_arclength, _point_distance, _closest_point_on_segment, _interpolate_point
import Proj

# Geodesic implementations
function GeometryOps._arclength_to_point(method::Geodesic, linestring, target_point, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    cumulative_distance = 0.0
    closest_distance = Inf
    result_distance = 0.0
    
    if GI.npoint(linestring) < 2
        return 0.0
    end
    
    proj_geodesic = Proj.geod_geodesic(method.semimajor_axis, 1/method.inv_flattening)
    prev_point = GI.getpoint(linestring, 1)
    
    for i in 2:GI.npoint(linestring)
        curr_point = GI.getpoint(linestring, i)
        
        # Calculate geodesic distance between consecutive points
        segment_length = _point_distance(method, prev_point, curr_point, proj_geodesic)
        
        # Find closest point on this geodesic segment to target
        closest_point_on_segment, t = _closest_point_on_segment(method, prev_point, curr_point, target_point, proj_geodesic)
        distance_to_segment = _point_distance(method, target_point, closest_point_on_segment, proj_geodesic)
        
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

function GeometryOps._point_at_arclength(method::Geodesic, linestring, target_distance, T::Union{GI.LineStringTrait, GI.LinearRingTrait})
    if GI.npoint(linestring) < 2
        return GI.npoint(linestring) > 0 ? GI.getpoint(linestring, 1) : nothing
    end
    
    if target_distance <= 0
        return GI.getpoint(linestring, 1)
    end
    
    proj_geodesic = Proj.geod_geodesic(method.semimajor_axis, 1/method.inv_flattening)
    cumulative_distance = 0.0
    prev_point = GI.getpoint(linestring, 1)
    
    for i in 2:GI.npoint(linestring)
        curr_point = GI.getpoint(linestring, i)
        segment_length = _point_distance(method, prev_point, curr_point, proj_geodesic)
        
        if cumulative_distance + segment_length >= target_distance
            # Target distance is within this segment
            remaining_distance = target_distance - cumulative_distance
            t = remaining_distance / segment_length
            return _interpolate_point(method, prev_point, curr_point, t, proj_geodesic)
        end
        
        cumulative_distance += segment_length
        prev_point = curr_point
    end
    
    # Distance exceeds total length, return last point
    return GI.getpoint(linestring, GI.npoint(linestring))
end

# Helper functions for geodesic distance calculations
function GeometryOps._point_distance(::Geodesic, p1, p2, proj_geodesic)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    geod_line = Proj.geod_inverseline(proj_geodesic, y1, x1, y2, x2)
    return geod_line.s13  # Distance in meters
end

# Find the closest point on a geodesic segment to a target point
function GeometryOps._closest_point_on_segment(::Geodesic, p1, p2, target, proj_geodesic)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    tx, ty = GI.x(target), GI.y(target)
    
    # Create geodesic line from p1 to p2
    geod_line = Proj.geod_inverseline(proj_geodesic, y1, x1, y2, x2)
    segment_length = geod_line.s13
    
    if segment_length == 0
        # Degenerate segment
        return p1, 0.0
    end
    
    # Sample points along the geodesic and find the closest to target
    best_t = 0.0
    min_distance = Inf
    
    # Sample at regular intervals to find approximate closest point
    n_samples = 100
    for i in 0:n_samples
        t = i / n_samples
        distance_along = t * segment_length
        
        # Get point at this position along geodesic
        lat, lon, _ = Proj.geod_position(geod_line, distance_along)
        sample_point = (lon, lat)
        
        # Calculate distance from sample point to target
        dist_to_target = _point_distance(Geodesic(), target, sample_point, proj_geodesic)
        
        if dist_to_target < min_distance
            min_distance = dist_to_target
            best_t = t
        end
    end
    
    # Get the closest point coordinates
    distance_along = best_t * segment_length
    lat, lon, _ = Proj.geod_position(geod_line, distance_along)
    closest_point = (lon, lat)
    
    return closest_point, best_t
end

# Interpolate between two points along a geodesic
function GeometryOps._interpolate_point(::Geodesic, p1, p2, t, proj_geodesic)
    x1, y1 = GI.x(p1), GI.y(p1)
    x2, y2 = GI.x(p2), GI.y(p2)
    
    # Create geodesic line from p1 to p2
    geod_line = Proj.geod_inverseline(proj_geodesic, y1, x1, y2, x2)
    distance_along = t * geod_line.s13
    
    # Get point at this position along geodesic
    lat, lon, _ = Proj.geod_position(geod_line, distance_along)
    
    return (lon, lat)
end