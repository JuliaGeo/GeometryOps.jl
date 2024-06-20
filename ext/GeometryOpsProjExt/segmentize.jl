# This holds the `segmentize` geodesic functionality.

import GeometryOps: GeodesicSegments, _fill_linear_kernel!, SVPoint_2D
import Proj

function GeometryOps.GeodesicSegments(; max_distance, equatorial_radius::Real=6378137, flattening::Real=1/298.257223563, geodesic::Proj.geod_geodesic = Proj.geod_geodesic(equatorial_radius, flattening))
    return GeometryOps.GeodesicSegments{Proj.geod_geodesic}(geodesic, max_distance)
end


function GeometryOps._fill_linear_kernel!(::Type{T}, method::GeodesicSegments{Proj.geod_geodesic}, new_coords::Vector, x1, y1, x2, y2) where T
    geod_line = Proj.geod_inverseline(method.geodesic, y1, x1, y2, x2)
    # This is the distance in meters computed between the two points.
    # It's `s13` because `geod_inverseline` sets point 3 to the second input point.
    distance = geod_line.s13 
    if distance > method.max_distance
        n_segments = ceil(Int, distance / method.max_distance)
        for i in 1:(n_segments - 1)
            y, x, _ = Proj.geod_position(geod_line, i / n_segments * distance)
            push!(new_coords, GeometryOps.SVPoint_2D((x, y), T))
        end
    end
    # End the line with the original coordinate,
    # to avoid any multiplication errors.
    push!(new_coords, GeometryOps.SVPoint_2D((x2, y2), T))
    return nothing
end