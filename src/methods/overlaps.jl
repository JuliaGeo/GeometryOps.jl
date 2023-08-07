# # Overlap checks

export overlaps

# This code checks whether geometries overlap with each other.  

# It does not compute the overlap or intersection geometry.

"""
    overlaps(geom1, geom2)::Bool

Compare two Geometries of the same dimension and return true if their intersection set results in a geometry
different from both but of the same dimension. It applies to Polygon/Polygon, LineString/LineString,
Multipoint/Multipoint, MultiLineString/MultiLineString and MultiPolygon/MultiPolygon.

## Examples
```jldoctest
julia> poly1 = Polygon([[[0,0],[0,5],[5,5],[5,0],[0,0]]])
Polygon(Array{Array{Float64,1},1}[[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]])

julia> poly2 = Polygon([[[1,1],[1,6],[6,6],[6,1],[1,1]]])
Polygon(Array{Array{Float64,1},1}[[[1.0, 1.0], [1.0, 6.0], [6.0, 6.0], [6.0, 1.0], [1.0, 1.0]]])

julia> overlap(poly1, poly2)
true
```
"""
overlaps(g1, g2)::Bool = overlaps(trait(g1), g1, trait(g2), g2)::Bool
overlaps(t1::FeatureTrait, g1, t2, g2)::Bool = overlaps(GI.geometry(g1), g2)
overlaps(t1, g1, t2::FeatureTrait, g2)::Bool = overlaps(g1, geometry(g2))
overlaps(t1::FeatureTrait, g1, t2::FeatureTrait, g2)::Bool = overlaps(geometry(g1), geometry(g2))
overlaps(::PolygonTrait, mp, ::MultiPolygonTrait, p)::Bool = overlaps(p, mp)
function overlaps(::MultiPointTrait, g1, ::MultiPointTrait, g2)::Bool 
    for p1 in GI.getpoint(g1)
        for p2 in GI.getpoint(g2)
            equals(p1, p2) && return true
        end
    end
end
function overlaps(::PolygonTrait, g1, ::PolygonTrait, g2)::Bool
    return line_intersects(g1, g2)
end
function overlaps(t1::MultiPolygonTrait, mp, t2::PolygonTrait, p1)::Bool
    for p2 in GI.getgeom(mp)
        overlaps(p1, thp2) && return true
    end
end
function overlaps(::MultiPolygonTrait, g1, ::MultiPolygonTrait, g2)::Bool
    for p1 in GI.getgeom(g1)
        overlaps(PolygonTrait(), mp, PolygonTrait(), p1) && return true
    end
end
