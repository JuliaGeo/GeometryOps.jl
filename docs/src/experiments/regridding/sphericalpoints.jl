import GeoInterface as GI
import GeometryOps as GO
import LinearAlgebra
import LinearAlgebra: dot, cross

## Get the area of a LinearRing with coordinates in radians
struct SphericalPoint{T <: Real}
	data::NTuple{3, T}
end
SphericalPoint(x, y, z) = SphericalPoint((x, y, z))

# define the 4 basic mathematical operators elementwise on the data tuple
Base.:+(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .+ q.data)
Base.:-(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .- q.data)
Base.:*(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .* q.data)
Base.:/(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data ./ q.data)
# Define sum on a SphericalPoint to sum across its data
Base.sum(p::SphericalPoint) = sum(p.data)

# define dot and cross products
LinearAlgebra.dot(p::SphericalPoint, q::SphericalPoint) = sum(p * q)
function LinearAlgebra.cross(a::SphericalPoint, b::SphericalPoint)
	a1, a2, a3 = a.data
    b1, b2, b3 = b.data
	SphericalPoint((a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1))
end

# Using Eriksson's formula for the area of spherical triangles: https://www.jstor.org/stable/2691141
# This melts down only when two points are antipodal, which in our case will not happen.
function _unit_spherical_triangle_area(a, b, c)
    #t = abs(dot(a, cross(b, c)))
    #t /= 1 + dot(b,c) + dot(c, a) + dot(a, b)
    t = abs(dot(a, (cross(b - a, c - a))) / dot(b + a, c + a))
    2*atan(t)
end

_lonlat_to_sphericalpoint(p) = _lonlat_to_sphericalpoint(GI.x(p), GI.y(p))
function _lonlat_to_sphericalpoint(lon, lat)
    lonsin, loncos = sincosd(lon)
    latsin, latcos = sincosd(lat)
    x = latcos * loncos
    y = latcos * lonsin
    z = latsin
    return SphericalPoint(x,y,z)
end



# Extend area to spherical

# TODO: make this the other way around, but that can wait.
GO.area(m::GO.Planar, geoms) = GO.area(geoms)

function GO.area(m::GO.Spherical, geoms)
    return GO.applyreduce(+, GI.PolygonTrait(), geoms) do poly
        GO.area(m, GI.PolygonTrait(), poly)
    end * ((-m.radius^2)/ 2) # do this after the sum, to increase accuracy and minimize calculations.
end

function GO.area(m::GO.Spherical, ::GI.PolygonTrait, poly)
    area = abs(_ring_area(m, GI.getexterior(poly)))
    for interior in GI.gethole(poly)
        area -= abs(_ring_area(m, interior))
    end
    return area
end

function _ring_area(m::GO.Spherical, ring)
    # convert ring to a sequence of SphericalPoints
    points = _lonlat_to_sphericalpoint.(GI.getpoint(ring))[1:end-1] # deliberately drop the closing point
    p1, p2, p3 = points[1], points[2], points[3]
    # For a spherical polygon, we can compute the area by splitting it into triangles
    # and summing their areas. We use the first point as a common vertex for all triangles.
    area = 0.0
    # Sum areas of triangles formed by first point and consecutive pairs of points
    np = length(points)
    for i in 1:np
        p1, p2, p3 = p2, p3, points[i]
        area += _unit_spherical_triangle_area(p1, p2, p3)
    end
    return area
end

function _ring_area(m::GO.Spherical, ring)
    # Convert ring points to spherical coordinates
    points = GO.tuples(ring).geom
    
    # Remove last point if it's the same as first (closed ring)
    if points[end] == points[1]
        points = points[1:end-1]
    end
    
    n = length(points)
    if n < 3
        return 0.0
    end

    area = 0.0

    # Use L'Huilier's formula to sum the areas of spherical triangles
    # formed by first point and consecutive pairs of points
    for i in 1:n
        p1, p2, p3 = points[mod1(i-1, n)], points[mod1(i, n)], points[mod1(i+1, n)]
        area += sind(GI.y(p2)) * (GI.x(p3) - GI.x(p1))
    end
    
    return area
end




# Test the area calculation
p1 = GI.Polygon([GI.LinearRing(Point2f[(0, 0), (1, 0), (0, 1), (0, 0)] .- (Point2f(0.5, 0.5),))])

GO.area(GO.Spherical(), p1)