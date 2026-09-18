import GeometryOps: _distance_point, _distance_segment

function GeometryOps._distance_auto(m::Geodesic, geom1, geom2, ::Type{T}, (longitude_scale, latitude_scale); kwargs...) where T
    longitude_scale == latitude_scale == 1 && return GeometryOps.distance(m, geom1, geom2, T; kwargs...)
    # Distance is a scalar, so the inputs can be converted to degrees without converting any result back.
    to_degrees(p) = (longitude_scale * p[1], latitude_scale * p[2])
    return GeometryOps.distance(m, GeometryOps.transform(to_degrees, geom1), GeometryOps.transform(to_degrees, geom2), T; kwargs...)
end

_proj_geodesic(m::Geodesic) = Proj.geod_geodesic(m.semimajor_axis, _flattening(m.inv_flattening))

function _distance_point(m::Geodesic, ::Type{T}, p1, p2) where T
    distance, _, _ = Proj.geod_inverse(_proj_geodesic(m), GI.y(p1), GI.x(p1), GI.y(p2), GI.x(p2))
    return T(distance)
end

#=
The distance from `p0` to a point at arc length `s` along the geodesic from `p1`
to `p2` is stationary where the geodesic towards `p0` is perpendicular to the
segment.  The cosine of the angle between the two azimuths is the negated
derivative of that distance, so its sign says whether `p0` is ahead of or
behind the point at `s`.

If `p0` is ahead of `p1` and behind `p2`, the closest point is the root of that
cosine, bracketed by the segment and found with the Illinois variant of regula
falsi.  Otherwise the closest point is an endpoint.
=#
function _distance_segment(m::Geodesic, ::Type{T}, p0, p1, p2) where T
    proj_geodesic = _proj_geodesic(m)
    lat0, lon0 = GI.y(p0), GI.x(p0)
    line = Proj.geod_inverseline(proj_geodesic, GI.y(p1), GI.x(p1), GI.y(p2), GI.x(p2))
    segment_length = line.s13

    function distance_and_direction(s)
        lat, lon, line_azimuth = Proj.geod_position(line, s)
        distance, azimuth, _ = Proj.geod_inverse(proj_geodesic, lat, lon, lat0, lon0)
        return distance, cosd(azimuth - line_azimuth)
    end

    d_lo, f_lo = distance_and_direction(zero(segment_length))
    d_hi, f_hi = distance_and_direction(segment_length)
    min_distance = min(d_lo, d_hi)
    (iszero(min_distance) || f_lo <= 0 || f_hi >= 0) && return T(min_distance)

    lo, hi = zero(segment_length), segment_length
    side = 0
    for _ in 1:100
        s = (lo * f_hi - hi * f_lo) / (f_hi - f_lo)
        d, f = distance_and_direction(s)
        min_distance = min(min_distance, d)
        (abs(f) <= 1e-12 || hi - lo <= 1e-9) && break
        if f > 0
            lo, f_lo = s, f
            side == 1 && (f_hi /= 2)
            side = 1
        else
            hi, f_hi = s, f
            side == -1 && (f_lo /= 2)
            side = -1
        end
    end
    return T(min_distance)
end
