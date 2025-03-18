#=
# Antimeridian Cutting
=#
export CutAtAntimeridianAndPoles # TODO: too wordy?

#=
This correction cuts the geometry at the antimeridian and the poles, and returns a fixed geometry.

The implementation is translated from the implementation in https://github.com/gadomski/antimeridian, 
which is a Python package.  Several ports of that algorithm exist in e.g. Go, Rust, etc.

At some point we will have to go in and clean up the implementation, remove all the hardcoded code,
and make it more efficient by using raw geointerface, and not allocating so much (perhaps, by passing allocs around).

But for right now it works, that's all I need.
=#

"""
    CutAtAntimeridianAndPoles(; kwargs...) <: GeometryCorrection

This correction cuts the geometry at the antimeridian and the poles, and returns a fixed geometry.
"""
Base.@kwdef struct CutAtAntimeridianAndPoles <: GeometryCorrection
    "The value at the north pole, in your angular units"
    northpole::Float64 = 90.0 # TODO not used!!!
    "The value at the south pole, in your angular units"
    southpole::Float64 = -90.0 # TODO not used!!!
    "The value at the left edge of the antimeridian, in your angular units"
    left::Float64 = -180.0 # TODO not used!!!
    "The value at the right edge of the antimeridian, in your angular units"
    right::Float64 = 180.0 # TODO not used!!!
    "The period of the cyclic / cylindrical coordinate system's x value, usually computed automatically so you don't have to provide it."
    period::Float64 = right - left # TODO not used!!!
    "If the polygon is known to enclose the north pole, set this to true"
    force_north_pole::Bool=false # TODO not used!!!
    "If the polygon is known to enclose the south pole, set this to true"
    force_south_pole::Bool=false # TODO not used!!!
    "If true, use the great circle method to find the antimeridian crossing, otherwise use the flat projection method.  There is no reason to set this to be off."
    great_circle = true
end

function Base.show(io::IO, cut::CutAtAntimeridianAndPoles)
    print(io, "CutAtAntimeridianAndPoles(left=$(cut.left), right=$(cut.right))")
end
Base.show(io::IO, ::MIME"text/plain", cut::CutAtAntimeridianAndPoles) = Base.show(io, cut)

application_level(::CutAtAntimeridianAndPoles) = TraitTarget(GI.PolygonTrait(), GI.LineStringTrait(), GI.MultiLineStringTrait(), GI.MultiPolygonTrait())

function (c::CutAtAntimeridianAndPoles)(trait::GI.AbstractTrait, geom)
    return _AntimeridianHelpers.cut_at_antimeridian(trait, geom)
end

module _AntimeridianHelpers

import GeoInterface as GI
import ..GeometryOps as GO

# Custom cross product for 3D tuples
function _cross(a::Tuple{Float64,Float64,Float64}, b::Tuple{Float64,Float64,Float64})
    return (
        a[2]*b[3] - a[3]*b[2],
        a[3]*b[1] - a[1]*b[3],
        a[1]*b[2] - a[2]*b[1]
    )
end

# Convert spherical degrees to cartesian coordinates
function spherical_degrees_to_cartesian(point::Tuple{Float64,Float64})::Tuple{Float64,Float64,Float64}
    lon, lat = point
    slon, clon = sincosd(lon)
    slat, clat = sincosd(lat)
    return (
        clon * clat,
        slon * clat,
        slat
    )
end

# Calculate crossing latitude using great circle method
function crossing_latitude_great_circle(start::Tuple{Float64,Float64}, endpoint::Tuple{Float64,Float64})::Float64
    # Convert points to 3D vectors
    p1 = spherical_degrees_to_cartesian(start)
    p2 = spherical_degrees_to_cartesian(endpoint)
    
    # Cross product defines plane through both points
    n1 = _cross(p1, p2)
    
    # Unit vector -Y defines meridian plane
    n2 = (0.0, -1.0, 0.0)
    
    # Intersection of planes defined by cross product
    intersection = _cross(n1, n2)
    norm = sqrt(sum(intersection .^ 2))
    intersection = intersection ./ norm
    
    # Convert back to spherical coordinates (just need latitude)
    round(asind(intersection[3]), digits=7)
end

# Calculate crossing latitude using flat projection method
function crossing_latitude_flat(start::Tuple{Float64,Float64}, endpoint::Tuple{Float64,Float64})::Float64
    lat_delta = endpoint[2] - start[2]
    if endpoint[1] > 0
        round(
            start[2] + (180.0 - start[1]) * lat_delta / (endpoint[1] + 360.0 - start[1]),
            digits=7
        )
    else
        round(
            start[2] + (start[1] + 180.0) * lat_delta / (start[1] + 360.0 - endpoint[1]),
            digits=7
        )
    end
end

# Main crossing latitude calculation function
function crossing_latitude(start::Tuple{Float64,Float64}, endpoint::Tuple{Float64,Float64}, great_circle::Bool)::Float64
    abs(start[1]) == 180 && return start[2]
    abs(endpoint[1]) == 180 && return endpoint[2]
    
    return great_circle ? crossing_latitude_great_circle(start, endpoint) : crossing_latitude_flat(start, endpoint)
end

# Normalize coordinates to ensure longitudes are between -180 and 180
function normalize_coords(coords::Vector{Tuple{Float64,Float64}})::Vector{Tuple{Float64,Float64}}
    normalized = deepcopy(coords)
    all_on_antimeridian = true
    
    for i in eachindex(normalized)
        point = normalized[i]
        prev_point = normalized[mod1(i-1, length(normalized))]
        
        if isapprox(point[1], 180)
            if abs(point[2]) != 90 && isapprox(prev_point[1], -180)
                normalized[i] = (-180.0, point[2])
            else
                normalized[i] = (180.0, point[2])
            end
        elseif isapprox(point[1], -180)
            if abs(point[2]) != 90 && isapprox(prev_point[1], 180)
                normalized[i] = (180.0, point[2])
            else
                normalized[i] = (-180.0, point[2])
            end
        else
            normalized[i] = (((point[1] + 180) % 360) - 180, point[2])
            all_on_antimeridian = false
        end
    end
    
    return all_on_antimeridian ? coords : normalized
end

# Segment a ring of coordinates at antimeridian crossings
function segment_ring(coords::Vector{Tuple{Float64,Float64}}, great_circle::Bool)::Vector{Vector{Tuple{Float64,Float64}}}
    segments = Vector{Vector{Tuple{Float64,Float64}}}()
    current_segment = Tuple{Float64,Float64}[]
    
    for i in 1:length(coords)-1
        start = coords[i]
        endpoint = coords[i+1]
        push!(current_segment, start)
        
        # Check for antimeridian crossing
        if (endpoint[1] - start[1] > 180) && (endpoint[1] - start[1] != 360)  # left crossing
            lat = crossing_latitude(start, endpoint, great_circle)
            push!(current_segment, (-180.0, lat))
            push!(segments, current_segment)
            current_segment = [(180.0, lat)]
        elseif (start[1] - endpoint[1] > 180) && (start[1] - endpoint[1] != 360)  # right crossing
            lat = crossing_latitude(endpoint, start, great_circle)
            push!(current_segment, (180.0, lat))
            push!(segments, current_segment)
            current_segment = [(-180.0, lat)]
        end
    end
    
    # Handle last point and segment
    if isempty(segments)
        return Vector{Vector{Tuple{Float64,Float64}}}()  # No crossings
    elseif coords[end] == segments[1][1]
        # Join polygons
        segments[1] = vcat(current_segment, segments[1])
    else
        push!(current_segment, coords[end])
        push!(segments, current_segment)
    end
    
    return normalize_coords.(segments)
end

# Check if a segment is self-closing
function is_self_closing(segment::Vector{Tuple{Float64,Float64}})::Bool
    is_right = segment[end][1] == 180
    return segment[1][1] == segment[end][1] && (
        (is_right && segment[1][2] > segment[end][2]) ||
        (!is_right && segment[1][2] < segment[end][2])
    )
end

# Build polygons from segments
function build_polygons(segments::Vector{Vector{Tuple{Float64,Float64}}})::Vector{GI.Polygon}
    isempty(segments) && return GI.Polygon[]
    
    segment = pop!(segments)
    is_right = segment[end][1] == 180
    candidates = Tuple{Union{Nothing,Int},Float64}[]
    
    if is_self_closing(segment)
        push!(candidates, (nothing, segment[1][2]))
    end
    
    for (i, s) in enumerate(segments)
        if s[1][1] == segment[end][1]
            if (is_right && s[1][2] > segment[end][2] && 
                (!is_self_closing(s) || s[end][2] < segment[1][2])) ||
               (!is_right && s[1][2] < segment[end][2] && 
                (!is_self_closing(s) || s[end][2] > segment[1][2]))
                push!(candidates, (i, s[1][2]))
            end
        end
    end
    
    # Sort candidates
    sort!(candidates, by=x->x[2], rev=!is_right)
    
    index = isempty(candidates) ? nothing : candidates[1][1]
    
    if !isnothing(index)
        # Join segments and recurse
        segment = vcat(segment, segments[index])
        segments[index] = segment
        return build_polygons(segments)
    else
        # Handle self-joining segment
        polygons = build_polygons(segments)
        if !all(p == segment[1] for p in segment)
            push!(polygons, GI.Polygon([segment]))
        end
        return polygons
    end
end

# Main function to cut a polygon at the antimeridian
cut_at_antimeridian(x) = cut_at_antimeridian(GI.trait(x), x)

function cut_at_antimeridian(
    ::GI.PolygonTrait,
    polygon::T;
    force_north_pole::Bool=false,
    force_south_pole::Bool=false,
    fix_winding::Bool=true,
    great_circle::Bool=true
) where {T}
    # Get exterior ring
    exterior = GO.tuples(GI.getexterior(polygon)).geom
    exterior = normalize_coords(exterior)
    
    # Segment the exterior ring
    segments = segment_ring(exterior, great_circle)
    
    if isempty(segments)
        # No antimeridian crossing
        if fix_winding && GO.isclockwise(GI.LinearRing(exterior))
            coord_vecs = reverse.(getproperty.(GO.tuples.(GI.getring(polygon)), :geom))
            return GI.Polygon(normalize_coords.(coord_vecs))
        end
        return polygon
    end
    
    # Handle holes
    holes = Vector{Vector{Tuple{Float64,Float64}}}()
    for hole_idx in 1:GI.nhole(polygon)
        hole = GO.tuples(GI.gethole(polygon, hole_idx)).geom
        hole_segments = segment_ring(hole, great_circle)
        
        if !isempty(hole_segments)
            if fix_winding
                unwrapped = [(x % 360, y) for (x, y) in hole]
                if !GO.isclockwise(GI.LineString(unwrapped))
                    hole_segments = segment_ring(reverse(hole), great_circle)
                end
            end
            append!(segments, hole_segments)
        else
            push!(holes, hole)
        end
    end
    
    # Build final polygons
    result_polygons = build_polygons(segments)
    
    # Add holes to appropriate polygons
    for hole in holes
        for (i, poly) in enumerate(result_polygons)
            if GO.contains(poly, GI.Point(hole[1]))
                rings = GI.getring(poly)
                push!(rings, hole)
                result_polygons[i] = GI.Polygon(rings)
                break
            end
        end
    end
    
    return length(result_polygons) == 1 ? result_polygons[1] : GI.MultiPolygon(result_polygons)
end

function cut_at_antimeridian(::GI.AbstractCurveTrait, line::T; great_circle::Bool=true) where {T}
    coords = GO.tuples(line).geom
    segments = segment_ring(coords, great_circle)
    
    if isempty(segments)
        return line
    else
        return GI.MultiLineString(segments)
    end
end


function cut_at_antimeridian(::GI.MultiPolygonTrait, x; kwargs...)
    results = GI.Polygon[]
    for poly in GI.getgeom(x)
        result = cut_at_antimeridian(GI.PolygonTrait(), poly; kwargs...)
        if result isa GI.Polygon
            push!(results, result)
        elseif result isa GI.MultiPolygon
            append!(results, result.geom)
        end
    end
    return GI.MultiPolygon(results)
end

function cut_at_antimeridian(::GI.MultiLineStringTrait, multiline::T; great_circle::Bool=true) where {T}
    linestrings = Vector{Vector{Tuple{Float64,Float64}}}()
    
    for line in GI.getgeom(multiline)
        fixed = cut_at_antimeridian(GI.LineStringTrait(), line; great_circle)
        if fixed isa GI.LineString
            push!(linestrings, GO.tuples(fixed).geom)
        else
            append!(linestrings, GO.tuples.(GI.getgeom(fixed)).geom)
        end
    end
    
    return GI.MultiLineString(linestrings)
end

end