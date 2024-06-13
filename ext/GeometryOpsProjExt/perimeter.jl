
function GeometryOps.GeodesicDistance(; equatorial_radius::Real=6378137, flattening::Real=1/298.257223563, geodesic::Proj.geod_geodesic = Proj.geod_geodesic(equatorial_radius, flattening))
    GeometryOps.GeodesicDistance{Proj.geod_geodesic}(geodesic)
end

function GeometryOps.point_distance(alg::GeometryOps.GeodesicDistance, p1, p2, ::Type{T}) where T <: Number
    lon1 = Base.convert(Float64, GI.x(p1))
    lat1 = Base.convert(Float64, GI.y(p1))
    lon2 = Base.convert(Float64, GI.x(p2))
    lat2 = Base.convert(Float64, GI.y(p2))

    dist, _azi1, _azi2 = Proj.geod_inverse(alg.geodesic, lon1, lat1, lon2, lat2)
    return T(dist)
end