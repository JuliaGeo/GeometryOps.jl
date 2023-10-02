# # Centroid

export centroid, centroid_and_length, centroid_and_area

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

## Implementation

This is the GeoInterface-compatible implementation.

First, we implement a wrapper method that dispatches to the correct
implementation based on the geometry trait. This is also used in the
implementation, since it's a lot less work! 

Note that if you call centroid on a LineString or LinearRing, the
centroid_and_length function will be called due to the weighting scheme
described above, while centroid_and_area is called for polygons and
multipolygons. However, centroid_and_area can still be called on a
LineString or LinearRing when they are closed, for example as the interior hole
of a polygon.

The helper functions centroid_and_length and centroid_and_area are made
availible just in case the user also needs the area or length to decrease
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
centroid(trait, geom) = centroid_and_area(trait, geom)[1]

"""
    centroid_and_length(geom)::(GI.Point, ::Real)

Returns the centroid and length of a given line/ring. Note this is only valid
for line strings and linear rings.
"""
centroid_and_length(geom) = centroid_and_length(GI.trait(geom), geom)

"""
    centroid_and_area(
        ::Union{GI.LineStringTrait, GI.LinearRingTrait}, 
        geom,
    )::(GI.Point, ::Real)

Returns the centroid and area of a given geom.
"""
centroid_and_area(geom) = centroid_and_area(GI.trait(geom), geom)

"""
    centroid_and_length(geom)::(GI.Point, ::Real)

Returns the centroid and length of a given line/ring. Note this is only valid
for line strings and linear rings.
"""
function centroid_and_length(
    ::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom,
)
    T = typeof(GI.x(GI.getpoint(geom, 1)))
    # Initialize starting values
    xcentroid = T(0)
    ycentroid = T(0)
    length = T(0)
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
    return (xcentroid, ycentroid), length
end

"""
    centroid_and_area(
        ::Union{GI.LineStringTrait, GI.LinearRingTrait},
        geom,
    )::(GI.Point, ::Real)

Returns the centroid and area of a given a line string or a linear ring.
Note that this is only valid if the line segment or linear ring is closed. 
"""
function centroid_and_area(
    ::Union{GI.LineStringTrait, GI.LinearRingTrait},
    geom,
)
    T = typeof(GI.x(GI.getpoint(geom, 1)))
    # Check that the geometry is closed
    @assert(
        GI.getpoint(geom, 1) == GI.getpoint(geom, GI.ngeom(geom)),
        "centroid_and_area should only be used with closed geometries"
    )
    # Initialize starting values
    xcentroid = T(0)
    ycentroid = T(0)
    area = T(0)
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
    return (xcentroid, ycentroid), abs(area)
end

"""
    centroid_and_area(::GI.PolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and area of a given polygon.
"""
function centroid_and_area(::GI.PolygonTrait, geom)
    # Exterior ring's centroid and area
    (xcentroid, ycentroid), area = centroid_and_area(GI.getexterior(geom))
    # Weight exterior centroid by area
    xcentroid *= area
    ycentroid *= area
    # Loop over any holes within the polygon
    for hole in GI.gethole(geom)
        # Hole polygon's centroid and area
        (xinterior, yinterior), interior_area = centroid_and_area(hole)
        # Accumulate the area component into `area`
        area -= interior_area
        # Weighted average of centroid components
        xcentroid -= xinterior * interior_area
        ycentroid -= yinterior * interior_area
    end
    xcentroid /= area
    ycentroid /= area
    return (xcentroid, ycentroid), area
end

"""
    centroid_and_area(::GI.MultiPolygonTrait, geom)::(GI.Point, ::Real)

Returns the centroid and area of a given multipolygon.
"""
function centroid_and_area(::GI.MultiPolygonTrait, geom)
    # First polygon's centroid and area
    (xcentroid, ycentroid), area = centroid_and_area(GI.getpolygon(geom, 1))
    # Weight first polygon's centroid by area
    xcentroid *= area
    ycentroid *= area
    # Loop over any polygons within the multipolygon
    for i in 2:GI.ngeom(geom) #poly in GI.getpolygon(geom)
        # Polygon centroid and area
        (xpoly, ypoly), poly_area = centroid_and_area(GI.getpolygon(geom, i))
        # Accumulate the area component into `area`
        area += poly_area
        # Weighted average of centroid components
        xcentroid += xpoly * poly_area
        ycentroid += ypoly * poly_area
    end
    xcentroid /= area
    ycentroid /= area
    return (xcentroid, ycentroid), area
end