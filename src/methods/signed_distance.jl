export signed_distance


Base.@propagate_inbounds euclid_distance(p1, p2) = sqrt((GeoInterface.x(p2)-GeoInterface.x(p1))^2 + (GeoInterface.y(p2)-GeoInterface.y(p1))^2)
euclid_distance(x1, y1, x2, y2) = sqrt((x2-x1)^2 + (y2-y1)^2)



" Distance from p0 to the line segment formed by p1 and p2.  Implementation from Turf.jl."
function _distance(p0, p1, p2)
    x0, y0 = GeoInterface.x(p0), GeoInterface.y(p0)
    x1, y1 = GeoInterface.x(p1), GeoInterface.y(p1)
    x2, y2 = GeoInterface.x(p2), GeoInterface.y(p2)

    if x1 < x2
        xfirst, yfirst = x1, y1
        xlast, ylast = x2, y2
    else
        xfirst, yfirst = x2, y2
        xlast, ylast = x1, y1
    end
    
    v = (xlast - xfirst, ylast - yfirst)
    w = (x0 - xfirst, y0 - yfirst)

    c1 = sum(w .* v)
    if c1 <= 0
        return euclid_distance(x0, y0, xfirst, yfirst)
    end

    c2 = sum(v .* v)

    if c2 <= c1
        return euclid_distance(x0, y0, xlast, ylast)
    end

    b2 = c1 / c2

    return euclid_distance(x0, y0, xfirst + (b2 * v[1]), yfirst + (b2 * v[2]))
end


function _distance(linestring, xy)
    mindist = typemax(Float64)
    N = GeoInterface.npoint(linestring)
    @assert N ≥ 3
    p1 = GeoInterface.getpoint(linestring, 1)
    p2 = p1
    
    for point_ind in 2:N
        p2 = GeoInterface.getpoint(linestring, point_ind)
        newdist = _distance(xy, p1, p2) 
        if newdist < mindist
            mindist = newdist
        end
        p1 = p2
    end

    return mindist
end

function signed_distance(::GeoInterface.PolygonTrait, poly, x, y)

    xy = (x, y)
    mindist = _distance(GeoInterface.getexterior(poly), xy)

    @inbounds for hole in GeoInterface.gethole(poly)
        newdist = _distance(hole, xy)
        if newdist < mindist
            mindist = newdist
        end
    end

    if GeoInterface.contains(poly, GeoInterface.convert(Base.parentmodule(typeof(poly)), (x, y)))
        return mindist
    else
        return -mindist
    end
end

function signed_distance(::GeoInterface.MultiPolygonTrait, multipoly, x, y)
    distances = signed_distance.(GeoInterface.getpolygon(multipoly), x, y)
    max_val, max_ind = findmax(distances)
    return max_val
end


"""
    signed_distance(geom, x::Real, y::Real)::Float64

Calculates the signed distance from the geometry `geom` to the point
defined by `(x, y)`.  Points within `geom` have a negative distance,
and points outside of `geom` have a positive distance.

If `geom` is a MultiPolygon, then this function returns the maximum distance 
to any of the polygons in `geom`.
"""
signed_distance(geom, x, y) = signed_distance(GeoInterface.geomtrait(geom), geom, x, y)
