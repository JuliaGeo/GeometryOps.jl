# # Centroid

export centroid

# These are all GeometryBasics.jl methods so far.
# They need to be converted to GeoInterface.

# The reason that there is a `centroid_and_signed_area` function, 
# is because in conputing the centroid, you end up computing the signed area.  

# In some computational geometry applications this may be a useful
# source of efficiency, so I added it here.

# However, it's totally fine to ignore this and not have this code path.  
# We simply need to decide on this.

function centroid(ls::GB.LineString{2, T}) where T
    centroid = Point{2, T}(0)
    total_area = T(0)
    if length(ls) == 1
        return sum(ls[1])/2
    end

    p0 = ls[1][1]

    for i in 1:(length(ls)-1)
        p1 = ls[i][2]
        p2 = ls[i+1][2]
        area = signed_area(p0, p1, p2)
        centroid = centroid .+ Point{2, T}((p0[1] + p1[1] + p2[1])/3, (p0[2] + p1[2] + p2[2])/3) * area
        total_area += area
    end
    return centroid ./ total_area
end

# a more optimized function, so we only calculate signed area once!
function centroid_and_signed_area(ls::GB.LineString{2, T}) where T
    centroid = Point{2, T}(0)
    total_area = T(0)
    if length(ls) == 1
        return sum(ls[1])/2
    end

    p0 = ls[1][1]

    for i in 1:(length(ls)-1)
        p1 = ls[i][2]
        p2 = ls[i+1][2]
        area = signed_area(p0, p1, p2)
        centroid = centroid .+ Point{2, T}((p0[1] + p1[1] + p2[1])/3, (p0[2] + p1[2] + p2[2])/3) * area
        total_area += area
    end
    return (centroid ./ total_area, total_area)
end

function centroid(poly::GB.Polygon{2, T}) where T
    exterior_centroid, exterior_area = centroid_and_signed_area(poly.exterior)

    total_area = exterior_area
    interior_numerator = Point{2, T}(0)
    for interior in poly.interiors
        interior_centroid, interior_area = centroid_and_signed_area(interior)
        total_area += interior_area
        interior_numerator += interior_centroid * interior_area
    end

    return (exterior_centroid * exterior_area - interior_numerator) / total_area

end

function centroid(multipoly::GB.MultiPolygon)
    centroids = centroid.(multipoly.polygons)
    areas = signed_area.(multipoly.polygons)
    areas ./= sum(areas)

    return sum(centroids .* areas) / sum(areas)
end


function centroid(rect::GB.Rect{N, T}) where {N, T}
    return Point{N, T}(rect.origin .- rect.widths ./ 2)
end

function centroid(sphere::GB.HyperSphere{N, T}) where {N, T}
    return sphere.center
end
