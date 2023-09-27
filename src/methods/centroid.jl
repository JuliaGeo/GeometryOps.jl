# # Centroid

export centroid

#=
## What is the centroid?

The centroid is the geometric center of a line string or area(s). Note that
the centroid does not need to be inside of a concave area.

To provide an example, consider this concave polygon in the shape of a 'C':
```@example cshape
using GeometryOps
using GeometryOps.GeometryBasics
using Makie
using CairoMakie

cshape = Polygon([
    Point(0,0), Point(0,3), Point(3,3), Point(3,2), Point(1,2),
    Point(1,1), Point(3,1), Point(3,0), Point(0,0),
])
f, a, p = poly(cshape; axis = (; aspect = DataAspect()))
```
Let's see what the centroid looks like:
```@example cshape
cent = centroid(cshape)
plot!(a, cent)
f
```
The points are ordered in a clockwise fashion, which means that the signed area
is positive.  If we reverse the order of the points, we get a negative area.

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait.

This is also used in the implementation, since it's a lot less work! 
=#
"""
    centroid(geom)::GI.Point

Returns the centroid of a given line segment, linear ring, polygon, or
mutlipolygon.
"""
centroid(geom) = centroid(GI.trait(geom), geom)

"""
    centroid(geom)::GI.Point

Returns the centroid of a given line segment.
"""
centroid(trait::GI.LineStringTrait, geom)  =
    centroid_and_length(trait, geom)[1]

"""
    centroid(geom)::GI.Point

Returns the centroid of a given linear ring, polygon, or
mutlipolygon.
"""
centroid(
    trait::Union{GI.LinearRingTrait, GI.PolygonTrait, GI.MultiPolygonTrait}, 
    geom
) = centroid_and_signed_area(trait, geom)[1]

"""
    centroid_and_length(geom)::GI.Point

Returns the centroid and length of a given geom. Note this is only valid for
line strings.
"""
centroid_and_length(geom) = centroid_and_length(GI.trait(geom), geom)

"""
    centroid_and_length(::GI.LineStringTrait, geom)::(GI.Point, ::Real)

Returns the centroid and length of a given line segment.
"""
function centroid_and_length(::GI.LineStringTrait, geom)
    FT = Float64
    # Initialize starting values
    centroid = GI.Point(FT(0), FT(0))
    length = FT(0)
    point₁ = GI.getpoint(geom, 1)
    # Loop over line segments of line string
    for point₂ in GI.getpoint(geom)
        # Calculate length of line segment
        length_component = sqrt(
            (GI.x(point₂) - GI.x(point₁))^2 +
            (GI.y(point₂) - GI.y(point₁))^2
        )
        # Accumulate the line segment length into `length`
        length += length_component
        # Weighted average of line segment centroids
        centroid = centroid .+ ((point₁ .+ point₂) .* (length_component / 2))
        # Advance the point buffer by 1 point to move to next line segment
        point₁ = point₂
    end
    return centroid ./= length, length
end

"""
    centroid_and_length(geom)::GI.Point

Returns the centroid and area of a given geom. Note this is only valid for
linear rings, polygons, and multipolygons.
"""
centroid_and_signed_area(geom) = centroid_and_signed_area(GI.trait(geom), geom)

"""
    centroid_and_signed_area(::GI.LinearRingTrait, geom)::(GI.Point, ::Real)

Returns the centroid and signed area of a given linear ring.
"""
function centroid_and_signed_area(::GI.LinearRingTrait, geom)
    FT = Float64
    # Initialize starting values
    centroid = GI.Point(FT(0), FT(0))
    area = FT(0)
    point₁ = GI.getpoint(geom, 1)
    # Loop over line segments of linear ring
    for point₂ in GI.getpoint(geom)
        area_component = GI.x(point₁) * GI.y(point₂) -
            GI.x(point₁) * GI.y(point₂)
        # Accumulate the area component into `area`
        area += area_component
        # Weighted average of centroid components
        centroid = centroid .+ (point₁ .+ point₂) .* area_component
        # Advance the point buffer by 1 point
        point₁ = point₂
    end
    area /= 2
    return centroid ./= 6area, area
end

"""
    centroid_and_signed_area(::GI.PolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and signed area of a given polygon.
"""
function centroid_and_signed_area(::GI.PolygonTrait, geom)
    # Exterior polygon centroid and area
    centroid, area = centroid_and_signed_area(
        GI.LinearRingTrait,
        GI.getexterior(geom),
    )
    # Weight exterior centroid by area
    centroid *= abs(signed_area)
    # Loop over any holes within the polygon
    for hole in GI.gethole(geom)
        # Hole polygon's centroid and area
        interior_centroid, interior_area = centroid_and_signed_area(hole)
        interior_area = abs(interior_area)
        # Accumulate the area component into `area`
        area -= interior_area
        # Weighted average of centroid components
        centroid = centroid .- (interior_centroid .* interior_area)
    end
    return centroid ./= abs(area), area
end

"""
    centroid_and_signed_area(::GI.MultiPolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and signed area of a given multipolygon.
"""
function centroid_and_signed_area(::GI.MultiPolygonTrait, geom)
    FT = Float64
    # Initialize starting values
    centroid = GI.Point(FT(0), FT(0))
    area = FT(0)
    # Loop over any polygons within the multipolygon
    for poly in GI.getpolygon(geom)
        # Polygon centroid and area
        poly_centroid, poly_area = centroid_and_signed_area(poly)
        # Accumulate the area component into `area`
        area += poly_area
        # Weighted average of centroid components
        centroid = centroid .+ (poly_centroid .* poly_area)
    end
    return centroid ./= area, area
end