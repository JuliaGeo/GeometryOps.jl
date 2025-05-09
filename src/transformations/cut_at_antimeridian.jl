# # Antimeridian cutting

export cut_at_antimeridian

using StaticArrays, LinearAlgebra
import GeoInterface as GI, GeoFormatTypes as GFT
import ..UnitSpherical

"""
    cut_at_antimeridian(geom; left_edge=-180.0, center_edge=0.0, right_edge=180.0, great_circle=true)

Cut a geometry along the antimeridian. Returns GeoInterface wrapper geometries.

This function handles geometries that cross the antimeridian by cutting them and 
returning appropriate parts. If a geometry does not cross the antimeridian, it is
returned unchanged.

## Parameters
- `geom`: The input geometry to be cut
- `left_edge`: The left edge of the antimeridian (default: -180.0)
- `center_edge`: The center edge of the Earth (default: 0.0)
- `right_edge`: The right edge of the antimeridian (default: 180.0)
- `great_circle`: Whether to compute meridian crossings on the sphere (default: true)

## Returns
- If the geometry crosses the antimeridian, a MultiPolygon or MultiLineString containing the cut parts
- If the geometry does not cross the antimeridian, the original geometry

## Example
```julia
using GeometryOps
using GeoInterface

# Create a polygon that crosses the antimeridian
poly = GeoInterface.Polygon([[
    (170.0, 40.0), (170.0, 50.0), (-170.0, 50.0), (-170.0, 40.0), (170.0, 40.0)
]])

# Cut the polygon at the antimeridian
cut_poly = cut_at_antimeridian(poly)
```
"""
function cut_at_antimeridian(geom; 
                           left_edge::Float64=-180.0, 
                           center_edge::Float64=0.0, 
                           right_edge::Float64=180.0,
                           great_circle::Bool=true)
    # Handle different geometry types
    return _cut_at_antimeridian(GI.trait(geom), geom, left_edge, center_edge, right_edge, great_circle)
end

# Default method for unsupported geometry types
function _cut_at_antimeridian(::GI.AbstractTrait, geom, left_edge, center_edge, right_edge, great_circle)
    throw(ArgumentError("Unsupported geometry type: $(typeof(geom))"))
end

# Method for point geometries (these can't cross the antimeridian, so return as-is)
function _cut_at_antimeridian(::GI.PointTrait, geom, left_edge, center_edge, right_edge, great_circle)
    return geom
end

# Method for LineString geometries
function _cut_at_antimeridian(::GI.LineStringTrait, geom, left_edge, center_edge, right_edge, great_circle)
    coords = collect(GI.getpoint(geom))
    
    # Convert to standard form (normalize longitudes)
    normalized_coords = normalize_coords(coords, left_edge, right_edge)
    
    # Segment the line string at the antimeridian
    segments = segment_coords(normalized_coords, left_edge, right_edge, great_circle)
    
    # If no segments were found, the line string doesn't cross the antimeridian
    if isempty(segments)
        return geom
    end
    
    # Create line strings from segments
    if length(segments) == 1
        return GI.LineString(segments[1])
    else
        return GI.MultiLineString(segments)
    end
end

# Method for Polygon geometries
function _cut_at_antimeridian(::GI.PolygonTrait, geom, left_edge, center_edge, right_edge, great_circle)
    # Get exterior ring coordinates
    exterior = GI.getexterior(geom)
    exterior_coords = collect(GI.getpoint(exterior))
    
    # Convert to standard form (normalize longitudes)
    normalized_coords = normalize_coords(exterior_coords, left_edge, right_edge)
    
    # Segment the polygon at the antimeridian
    segments = segment_coords(normalized_coords, left_edge, right_edge, great_circle)
    
    # If no segments were found, the polygon doesn't cross the antimeridian
    if isempty(segments)
        return geom
    end
    
    # Extend segments over poles if needed
    segments = extend_segments_over_poles(segments, left_edge, right_edge)
    
    # Build polygons from segments
    polygons = build_polygons_from_segments(segments)
    
    # Handle interior rings (holes)
    holes = collect(GI.gethole(geom))
    if !isempty(holes)
        for interior in holes
            interior_coords = collect(GI.getpoint(interior))
            interior_norm = normalize_coords(interior_coords, left_edge, right_edge)
            interior_segments = segment_coords(interior_norm, left_edge, right_edge, great_circle)
            
            # If the interior ring crosses the antimeridian, add its segments
            if !isempty(interior_segments)
                # Add interior segments to appropriate polygons
                add_interior_segments_to_polygons!(polygons, interior_segments)
            else
                # If the interior doesn't cross, add it as a hole to the correct polygon
                add_interior_to_polygons!(polygons, interior, interior_coords)
            end
        end
    end
    
    # Create result as a single polygon or multipolygon
    if length(polygons) == 1
        return polygons[1]
    else
        return GI.MultiPolygon(polygons)
    end
end

# Method for MultiPolygon geometries
function _cut_at_antimeridian(::GI.MultiPolygonTrait, geom, left_edge, center_edge, right_edge, great_circle)
    # Process each polygon in the multipolygon
    result_polygons = []
    
    for poly in GI.getgeom(geom)
        # Cut each polygon
        cut_poly = _cut_at_antimeridian(GI.PolygonTrait(), poly, left_edge, center_edge, right_edge, great_circle)
        
        # Add the result to our collection
        if GI.geomtrait(cut_poly) isa GI.PolygonTrait
            push!(result_polygons, cut_poly)
        elseif GI.geomtrait(cut_poly) isa GI.MultiPolygonTrait
            append!(result_polygons, GI.getgeom(cut_poly))
        end
    end
    
    # Return as a MultiPolygon
    return GI.MultiPolygon(result_polygons)
end

# Method for MultiLineString geometries
function _cut_at_antimeridian(::GI.MultiLineStringTrait, geom, left_edge, center_edge, right_edge, great_circle)
    # Process each line string in the multilinestring
    result_linestrings = []
    
    for line in GI.getgeom(geom)
        # Cut each line string
        cut_line = _cut_at_antimeridian(GI.LineStringTrait(), line, left_edge, center_edge, right_edge, great_circle)
        
        # Add the result to our collection
        if GI.geomtrait(cut_line) isa GI.LineStringTrait
            push!(result_linestrings, cut_line)
        elseif GI.geomtrait(cut_line) isa GI.MultiLineStringTrait
            append!(result_linestrings, GI.getgeom(cut_line))
        end
    end
    
    # Return as a MultiLineString
    return GI.MultiLineString(result_linestrings)
end

# Normalize coordinates to standard range
function normalize_coords(coords, left_edge, right_edge)
    normalized = similar(coords)
    span = right_edge - left_edge
    
    for i in eachindex(coords)
        x, y = coords[i]
        # Normalize longitude to the range [left_edge, right_edge]
        x_norm = x
        if x < left_edge || x > right_edge
            x_norm = ((x - left_edge) % span + span) % span + left_edge
        end
        normalized[i] = (x_norm, y)
    end
    
    return normalized
end

# Segment coordinates at the antimeridian
function segment_coords(coords, left_edge, right_edge, great_circle)
    segments = Vector{Vector{Tuple{Float64, Float64}}}()
    current_segment = Tuple{Float64, Float64}[]
    span = right_edge - left_edge
    
    # Check pairs of coordinates for antimeridian crossing
    for i in 1:length(coords)-1
        start, finish = coords[i], coords[i+1]
        push!(current_segment, start)
        
        # Check for crossing from left to right
        if finish[1] - start[1] > span/2
            # Calculate the latitude at crossing
            lat = crossing_latitude(start, finish, left_edge, right_edge, great_circle)
            push!(current_segment, (left_edge, lat))
            push!(segments, current_segment)
            current_segment = [(right_edge, lat)]
        
        # Check for crossing from right to left
        elseif start[1] - finish[1] > span/2
            # Calculate the latitude at crossing
            lat = crossing_latitude(start, finish, left_edge, right_edge, great_circle)
            push!(current_segment, (right_edge, lat))
            push!(segments, current_segment)
            current_segment = [(left_edge, lat)]
        end
    end
    
    # Add the last point to the current segment
    push!(current_segment, coords[end])
    
    # If we found segments, add the last segment
    if !isempty(segments)
        push!(segments, current_segment)
    end
    
    return segments
end

# Calculate the latitude where a line segment crosses the antimeridian
function crossing_latitude(start, finish, left_edge, right_edge, great_circle)
    # Handle cases where points are already on the antimeridian
    if abs(start[1] - left_edge) < eps() || abs(start[1] - right_edge) < eps()
        return start[2]
    elseif abs(finish[1] - left_edge) < eps() || abs(finish[1] - right_edge) < eps()
        return finish[2]
    end
    
    if great_circle
        return crossing_latitude_great_circle(start, finish, left_edge, right_edge)
    else
        return crossing_latitude_flat(start, finish, left_edge, right_edge)
    end
end

# Calculate crossing latitude using great circle path (spherical)
function crossing_latitude_great_circle(start, finish, left_edge, right_edge)
    # Convert degrees to 3D Cartesian coordinates on the unit sphere
    start_point = UnitSpherical.UnitSphereFromGeographic()(start)
    finish_point = UnitSpherical.UnitSphereFromGeographic()(finish)
    
    # The meridian plane is defined by (0, -1, 0) for the -180/180 meridian
    # We need to adjust this based on left_edge and right_edge
    meridian_adjustment = (left_edge + right_edge) / 2
    meridian_rad = deg2rad(meridian_adjustment)
    meridian_normal = SVector{3}(sin(meridian_rad), -cos(meridian_rad), 0.0)
    
    # The plane containing the two points and the origin
    plane_normal = LinearAlgebra.cross(SVector{3}(start_point.x, start_point.y, start_point.z), 
                                       SVector{3}(finish_point.x, finish_point.y, finish_point.z))
    
    # The intersection of both planes is defined by their cross product
    intersection = LinearAlgebra.cross(plane_normal, meridian_normal)
    intersection = intersection / LinearAlgebra.norm(intersection)
    
    # Convert back to geographic coordinates
    point_on_sphere = UnitSpherical.UnitSphericalPoint(intersection[1], intersection[2], intersection[3])
    lon_lat = UnitSpherical.GeographicFromUnitSphere()(point_on_sphere)
    
    # Return just the latitude
    return lon_lat[2]
end

# Calculate crossing latitude using flat (Cartesian) path
function crossing_latitude_flat(start, finish, left_edge, right_edge)
    lat_delta = finish[2] - start[2]
    lon_delta = (finish[1] - start[1] + 360) % 360
    
    if lon_delta > 180
        lon_delta = 360 - lon_delta
        lat_delta = -lat_delta
    end
    
    if finish[1] > start[1]
        # Crossing from left to right
        proportion = (right_edge - start[1]) / lon_delta
    else
        # Crossing from right to left
        proportion = (start[1] - left_edge) / lon_delta
    end
    
    return start[2] + proportion * lat_delta
end

# Extend segments over poles if needed
function extend_segments_over_poles(segments, left_edge, right_edge)
    # Implementation based on Python's extend_over_poles
    # This function identifies segments that should be connected over poles
    
    # If we have fewer than 2 segments, no need to extend
    if length(segments) < 2
        return segments
    end
    
    # Find segments that end at the same latitude at the antimeridian
    # These are candidates for connecting over poles
    extended_segments = copy(segments)
    
    # Group segments by their endpoint latitudes
    left_endpoints = Dict{Float64, Vector{Int}}()
    right_endpoints = Dict{Float64, Vector{Int}}()
    
    for (i, segment) in enumerate(segments)
        # Check if segment starts at left edge
        if abs(segment[1][1] - left_edge) < eps()
            lat = segment[1][2]
            if !haskey(left_endpoints, lat)
                left_endpoints[lat] = []
            end
            push!(left_endpoints[lat], i)
        end
        
        # Check if segment ends at right edge
        if abs(segment[end][1] - right_edge) < eps()
            lat = segment[end][2]
            if !haskey(right_endpoints, lat)
                right_endpoints[lat] = []
            end
            push!(right_endpoints[lat], i)
        end
    end
    
    # Connect segments over poles if they have matching endpoints
    for (lat, left_indices) in left_endpoints
        if haskey(right_endpoints, lat)
            for left_idx in left_indices
                for right_idx in right_endpoints[lat]
                    # Connect these segments
                    # For simplicity, we'll just merge them directly
                    # In a more sophisticated implementation, we might add points along the pole
                    left_segment = extended_segments[left_idx]
                    right_segment = extended_segments[right_idx]
                    
                    # Create a new segment that combines both
                    new_segment = vcat(right_segment, left_segment[2:end])
                    
                    # Replace the right segment with the combined one
                    extended_segments[right_idx] = new_segment
                    
                    # Mark the left segment for removal
                    extended_segments[left_idx] = Tuple{Float64, Float64}[]
                end
            end
        end
    end
    
    # Remove empty segments
    filter!(s -> !isempty(s), extended_segments)
    
    return extended_segments
end

# Build polygons from segments
function build_polygons_from_segments(segments)
    # Implementation based on Python's build_polygons
    # For each segment, create a polygon
    polygons = []
    
    for segment in segments
        # Ensure the segment is closed
        if segment[1] != segment[end]
            push!(segment, segment[1])
        end
        push!(polygons, GI.Polygon([segment]))
    end
    
    return polygons
end

# Add interior segments to appropriate polygons
function add_interior_segments_to_polygons!(polygons, interior_segments)
    # This would assign interior segments to the right polygons
    # For simplicity in this implementation, we'll just add them as separate polygons
    # with a "hole" flag
    
    for segment in interior_segments
        # Ensure the segment is closed
        if segment[1] != segment[end]
            push!(segment, segment[1])
        end
        
        # Find which polygon contains this interior segment
        for (i, poly) in enumerate(polygons)
            # For now, we'll just add it to the first polygon
            # In a more sophisticated implementation, we'd check containment
            exterior = GI.getexterior(poly)
            holes = collect(GI.gethole(poly))
            push!(holes, segment)
            
            # Replace the polygon with a new one that includes the hole
            polygons[i] = GI.Polygon([exterior, holes...])
            break
        end
    end
end

# Add an interior ring to the appropriate polygon
function add_interior_to_polygons!(polygons, interior, interior_coords)
    # This adds the interior as a hole to the polygon that contains it
    
    # Find which polygon contains this interior
    for (i, poly) in enumerate(polygons)
        # For simplicity, we'll just add it to the first polygon
        # In a more sophisticated implementation, we'd check containment
        exterior = GI.getexterior(poly)
        holes = collect(GI.gethole(poly))
        push!(holes, interior)
        
        # Replace the polygon with a new one that includes the hole
        polygons[i] = GI.Polygon([exterior, holes...])
        break
    end
end
