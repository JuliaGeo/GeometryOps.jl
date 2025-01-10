# This holds the `segmentize` geodesic functionality.

import GeometryOps: GeodesicSegments, _segmentize, _fill_linear_kernel!
import Proj

function GeometryOps.GeodesicSegments(; max_distance, equatorial_radius::Real=6378137, flattening::Real=1/298.257223563, geodesic::Proj.geod_geodesic = Proj.geod_geodesic(equatorial_radius, flattening))
    return GeometryOps.GeodesicSegments{Proj.geod_geodesic}(geodesic, max_distance)
end

function _normalize_longitudes(data)
    GeometryOps.transform(data) do p
        x, y = GI.x(p), GI.y(p)
        return GI.Point(x > 180 ? x - 360 : x, y)
    end
end

# This is the same method as in `transformations/segmentize.jl`,
# but it constructs a Proj geodesic line every time.
# Maybe this should be better...
function _segmentize(method::Geodesic, geom, ::Union{GI.LineStringTrait, GI.LinearRingTrait}; max_distance)
    # Construct the Proj geodesic reference.
    proj_geodesic = Proj.geod_geodesic(method.semimajor_axis #= same thing as equatorial radius =#, 1/method.inv_flattening)
    # Get and normalize the first coordinate, so we don't break with Proj's assumptions.
    first_coord = GI.getpoint(geom, 1)
    second_coord = GI.getpoint(geom, 2)
    x1, y1 = GI.x(first_coord), GI.y(first_coord)
    x2, y2 = GI.x(second_coord), GI.y(second_coord)
    geod_line = Proj.geod_inverseline(proj_geodesic, y1, x1, y2, x2)
    y1_norm, x1_norm, _ = Proj.geod_position(geod_line, 0.0)

    new_coords = NTuple{2, Float64}[]
    sizehint!(new_coords, GI.npoint(geom))

    push!(new_coords, (x1_norm, y1_norm))
    # Iterate over the rest of the coordinates.
    for coord in Iterators.drop(GI.getpoint(geom), 1)
        x2, y2 = GI.x(coord), GI.y(coord)
        _fill_linear_kernel!(method, new_coords, x1, y1, x2, y2; max_distance, proj_geodesic)
        x1, y1 = new_coords[end]
    end 
    return rebuild(geom, new_coords)
end

function GeometryOps._fill_linear_kernel!(method::Geodesic, new_coords::Vector, x1, y1, x2, y2; max_distance, proj_geodesic)
    geod_line = Proj.geod_inverseline(proj_geodesic, y1, x1, y2, x2)
    # This is the distance in meters computed between the two points.
    # It's `s13` because `geod_inverseline` sets point 3 to the second input point.
    distance = geod_line.s13 
    if distance > max_distance
        n_segments = ceil(Int, distance / max_distance)
        for i in 1:(n_segments - 1)
            y, x, _ = Proj.geod_position(geod_line, i / n_segments * distance)
            push!(new_coords, (x, y))
        end
    end
    # End the line with the original coordinate,
    # to avoid any multiplication errors.
    yf, xf, _ = Proj.geod_position_relative(geod_line, 1.0)
    push!(new_coords, (xf, yf))
    return nothing
end