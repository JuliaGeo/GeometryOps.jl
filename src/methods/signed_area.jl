# # Signed area
export signed_area

# This is the GeoInterface-compatible implementation.

# First, we implement a wrapper method that dispatches to the correct
# implementation based on the geometry trait.
# 
# This is also used in the implementation, since it's a lot less work! 
"""
    signed_area(geom)::Real

Returns the signed area of the geometry, based on winding order.
"""
signed_area(x) = signed_area(GI.trait(x), x)

# TODOS here:
# 1. This could conceivably be multithreaded.  How to indicate that it should be so?
# 2. What to do for corner cases (nan point, etc)?
function signed_area(::Union{LineStringTrait, LinearRingTrait}, geom)
    # Basically, we integrate the area under the line string, which gives us
    # the signed area.
    point₁ = GI.getpoint(geom, 1)
    point₂ = GI.getpoint(geom, 2)
    area = GI.x(point₁) * GI.y(point₂) - GI.y(point₁) * GI.x(point₂)
    for point in GI.getpoint(geom)
        # Advance the point buffers by 1 point
        point₁ = point₂
        point₂ = point
        # Accumulate the area into `area`
        area += GI.x(point₁) * GI.y(point₂) - GI.y(point₁) * GI.x(point₂)
    end
    area /= 2
    return area
end

# This subtracts the 
function signed_area(::PolygonTrait, geom)
    s_area = signed_area(GI.getexterior(geom))
    area = abs(s_area)
    for hole in GI.gethole(geom)
        area -= abs(signed_area(hole))
    end
    return area * sign(s_area)
end

signed_area(::MultiPolygonTrait, geom) = sum((signed_area(poly) for poly in GI.getpolygon(geom)))

# This should _theoretically_ work for anything, but I haven't actually tested yet!

# Below is the original GeometryBasics implementation:

# # ```julia
# function signed_area(a::Point{2, T}, b::Point{2, T}, c::Point{2, T}) where T
#     return ((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
# end
#
# function signed_area(points::AbstractVector{<: Point{2, T}}) where {T}
#     area = sum((points[i][1] * points[i+1][2] - points[i][2] * points[i+1][1] for i in 1:(length(points)-1))) / 2.0
# end
#
# function signed_area(ls::GeometryBasics.LineString)
#     # coords = GeometryBasics.decompose(Point2f, ls)
#     return sum((p1[1] * p2[2] - p1[2] * p2[1] for (p1, p2) in ls)) / 2.0#signed_area(coords)
# end
#
# function signed_area(poly::GeometryBasics.Polygon{2})
#     s_area = signed_area(poly.exterior)
#     area = abs(s_area)
#     for hole in poly.interiors
#         area -= abs(signed_area(hole))
#     end
#     return area * sign(s_area)
# end
#
# # WARNING: this may not do what you expect, since it's
# # sensitive to winding order.  Use GeoInterface.area instead.
# signed_area(mp::MultiPolygon) = sum(signed_area.(mp.polygons))
# ```