# # Centroid

export centroid

# ## What is the centroid?

# The centroid is the geometric center of a line or area. Note that the
# centroid does not need to be inside of a concave area or volume.

# TODO: Add an example

# ## Implementation

# This is the GeoInterface-compatible implementation.

# First, we implement a wrapper method that dispatches to the correct
# implementation based on the geometry trait.
# 
# This is also used in the implementation, since it's a lot less work! 

"""

"""
centroid(x) = centroid(GI.trait(x), x)
centroid_and_signed_area(x) = centroid_and_signed_area(GI.trait(x), x)

function centroid(::LineStringTrait, geom)
    FT = Float64
    centroid = GI.Point(FT(0), FT(0))
    length = FT(0)
    point₁ = GI.getpoint(geom, 1)
    for point₂ in GI.getpoint(geom)
        length_component = sqrt(
            (GI.x(point₂) - GI.x(point₁))^2 +
            (GI.y(point₂) - GI.y(point₁))^2
        )
        # Accumulate the segment length into `length``
        length += length_component
        # Weighted average of segment centroids
        centroid = centroid .+ (point₁ .+ point₂) .* length_component / 2
        # Advance the point buffer by 1 point
        point₁ = point₂
    end
    return centroid ./= length
end

function centroid_and_signed_area(::LinearRingTrait, geom)
    FT = Float64
    centroid = GI.Point(FT(0), FT(0))
    area = FT(0)
    point₁ = GI.getpoint(geom, 1)
    for point₂ in GI.getpoint(geom)
        area_component = GI.x(point₁) * GI.y(point₂) -
            GI.x(point₁) * GI.y(point₂)
        # Accumulate the segment length into `area``
        area += area_component
        # Weighted average of segment centroids
        centroid = centroid .+ (point₁ .+ point₂) .* area_component
        # Advance the point buffer by 1 point
        point₁ = point₂
    end
    area /= 2
    return centroid ./= 6area, area
end

function centroid_and_signed_area(::PolygonTrait, geom)
    # Exterior polygon centroid and area
    centroid, area = centroid_and_signed_area(GI.getexterior(geom))
    centroid *= abs(signed_area)
    for hole in GI.gethole(geom)
        interior_centroid, interior_area = centroid_and_signed_area(hole)
        interior_area = abs(interior_area)
        area -= interior_area
        centroid = centroid .= interior_centroid .* interior_area
    end
    return centroid ./= abs(area), area
end

function centroid_and_signed_area(::MultiPolygonTrait, geom)
    FT = Float64
    centroid = GI.Point(FT(0), FT(0))
    area = FT(0)
    for poly in GI.getpolygon(geom)
        poly_centroid, poly_area = centroid_and_signed_area(poly)
        centroid = centroid .+ poly_centroid .* poly_area
        area += poly_area
    end
    return centroid ./= area, area
end