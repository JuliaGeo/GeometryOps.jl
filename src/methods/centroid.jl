# # Centroid

export centroid

#=
## What is the centroid?

The centroid is the geometric center of a line string or area(s). Note that
the centroid does not need to be inside of a concave area.

Further note that by convention a line, or linear ring, is calculated by
weighting the line segments by their length, while polygons and multipolygon
centroids are calculated by weighting edge's by their 'area components'.

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
Let's see what the centroid looks like (plotted in red):
```@example cshape
cent = centroid(cshape)
scatter!(a, GI.x(cent), GI.y(cent), color = :red)
f
```
The points are ordered in a clockwise fashion, which means that the signed area
is positive.  If we reverse the order of the points, we get a negative area.

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait. This is also used in the
implementation, since it's a lot less work! 

Note that if you call centroid on a LineString or LinearRing, the
centroid_and_length function will be called due to the weighting scheme
described above, while centroid_and_signed_area is called for polygons and
multipolygons. However, centroid_and_signed_area can still be called on a
LineString or LinearRing when they are closed, for example as the interior hole
of a polygon.

The helper functions centroid_and_length and centroid_and_signed_area are made
availible just in case the user also needs the signed area or length to decrease
repeat computation.
=#
"""
    centroid(geom)::GI.Point

Returns the centroid of a given line segment, linear ring, polygon, or
mutlipolygon.
"""
centroid(geom) = centroid(GI.trait(geom), geom)

"""
    centroid(
        trait::Union{GI.LineStringTrait, GI.LinearRingTrait},
        geom,
    )

Returns the centroid of a line string or linear ring, which is calculated by
weighting line segments by their length by convention.
"""
centroid(
    trait::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom,
) = centroid_and_length(trait, geom)[1]

"""
    centroid(trait, geom)

Returns the centroid of a polygon or multipolygon, which is calculated by
weighting edges by their `area component` by convention.
"""
centroid(trait, geom) = centroid_and_signed_area(trait, geom)[1]

"""
    centroid_and_length(geom)::(GI.Point, ::Real)

Returns the centroid and length of a given line/ring. Note this is only valid
for line strings and linear rings.
"""
centroid_and_length(geom) = centroid_and_length(GI.trait(geom), geom)

"""
    centroid_and_signed_area(
        ::Union{GI.LineStringTrait, GI.LinearRingTrait}, 
        geom,
    )::(GI.Point, ::Real)

Returns the centroid and signed area of a given geom.
"""
centroid_and_signed_area(geom) = centroid_and_signed_area(GI.trait(geom), geom)

"""
    centroid_and_length(geom)::(GI.Point, ::Real)

Returns the centroid and length of a given line/ring. Note this is only valid
for line strings and linear rings.
"""
function centroid_and_length(
    ::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom,
)
    FT = Float64
    # Initialize starting values
    xcentroid = FT(0)
    ycentroid = FT(0)
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
        xcentroid += (GI.x(point₁) + GI.x(point₂)) * (length_component / 2)
        ycentroid += (GI.y(point₁) + GI.y(point₂)) * (length_component / 2)
        #centroid = centroid .+ ((point₁ .+ point₂) .* (length_component / 2))
        # Advance the point buffer by 1 point to move to next line segment
        point₁ = point₂
    end
    xcentroid /= length
    ycentroid /= length
    return GI.Point(xcentroid, ycentroid), length
end

"""
    centroid_and_signed_area(
        ::Union{GI.LineStringTrait, GI.LinearRingTrait},
        geom,
    )::(GI.Point, ::Real)

Returns the centroid and signed area of a given a line string or a linear ring.
Note that the area doesn't have much meaning as for a line string, and isn't
valid if the line segment isn't closed. 
"""
function centroid_and_signed_area(
    ::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom,
)
    FT = Float64
    # Initialize starting values
    xcentroid = FT(0)
    ycentroid = FT(0)
    area = FT(0)
    point₁ = GI.getpoint(geom, 1)
    # Loop over line segments of linear ring
    for point₂ in GI.getpoint(geom)
        area_component = GI.x(point₁) * GI.y(point₂) -
            GI.x(point₂) * GI.y(point₁)
        # Accumulate the area component into `area`
        area += area_component
        # Weighted average of centroid components
        xcentroid += (GI.x(point₁) + GI.x(point₂)) * area_component
        ycentroid += (GI.y(point₁) + GI.y(point₂)) * area_component
        # Advance the point buffer by 1 point
        point₁ = point₂
    end
    area /= 2
    xcentroid /= 6area
    ycentroid /= 6area
    return GI.Point(xcentroid, ycentroid), area
end

"""
    centroid_and_signed_area(::GI.PolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and signed area of a given polygon.
"""
function centroid_and_signed_area(::GI.PolygonTrait, geom)
    FT = Float64
    # Initialize starting values
    xcentroid = FT(0)
    ycentroid = FT(0)
    area = FT(0)
    # Exterior polygon centroid and area
    ext_centroid, ext_area = centroid_and_signed_area(GI.getexterior(geom))
    area += ext_area
    ext_area = abs(ext_area)
    # Weight exterior centroid by area
    xcentroid += GI.x(ext_centroid) * ext_area
    ycentroid += GI.y(ext_centroid) * ext_area
    # Loop over any holes within the polygon
    for hole in GI.gethole(geom)
        # Hole polygon's centroid and area
        interior_centroid, interior_area = centroid_and_signed_area(hole)
        interior_area = abs(interior_area)
        # Accumulate the area component into `area`
        area -= interior_area
        # Weighted average of centroid components
        xcentroid -= GI.x(interior_centroid) * interior_area
        ycentroid -= GI.y(interior_centroid) * interior_area
    end
    xcentroid /= abs(area)
    ycentroid /= abs(area)
    return GI.Point(xcentroid, ycentroid), area
end

"""
    centroid_and_signed_area(::GI.MultiPolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and signed area of a given multipolygon.
"""
function centroid_and_signed_area(::GI.MultiPolygonTrait, geom)
    FT = Float64
    # Initialize starting values
    xcentroid = FT(0)
    ycentroid = FT(0)
    area = FT(0)
    # Loop over any polygons within the multipolygon
    for poly in GI.getpolygon(geom)
        # Polygon centroid and area
        poly_centroid, poly_area = centroid_and_signed_area(poly)
        # Accumulate the area component into `area`
        area += poly_area
        # Weighted average of centroid components
        xcentroid += GI.x(poly_centroid) * poly_area
        ycentroid += GI.y(poly_centroid) * poly_area
    end
    xcentroid /= area
    ycentroid /= area
    return GI.Point(xcentroid, ycentroid), area
end