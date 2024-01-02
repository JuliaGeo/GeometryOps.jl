# # Crossing checks

"""
     crosses(geom1, geom2)::Bool

Return `true` if the intersection results in a geometry whose dimension is one less than
the maximum dimension of the two source geometries and the intersection set is interior to
both source geometries.

TODO: broken

## Examples 
```julia
import GeoInterface as GI, GeometryOps as GO

line1 = GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])
line2 = GI.LineString([(-2, 2), (4, 2)])

GO.crosses(line1, line2)
# output
true
```
"""
crosses(g1, g2)::Bool = crosses(trait(g1), g1, trait(g2), g2)::Bool
crosses(t1::FeatureTrait, g1, t2, g2)::Bool = crosses(GI.geometry(g1), g2)
crosses(t1, g1, t2::FeatureTrait, g2)::Bool = crosses(g1, geometry(g2))
crosses(::MultiPointTrait, g1, ::LineStringTrait, g2)::Bool = multipoint_crosses_line(g1, g2)
crosses(::MultiPointTrait, g1, ::PolygonTrait, g2)::Bool = multipoint_crosses_poly(g1, g2)
crosses(::LineStringTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_crosses_lines(g2, g1)
crosses(::LineStringTrait, g1, ::PolygonTrait, g2)::Bool = line_crosses_poly(g1, g2)
crosses(::LineStringTrait, g1, ::LineStringTrait, g2)::Bool = line_crosses_line(g1, g2)
crosses(::PolygonTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_crosses_poly(g2, g1)
crosses(::PolygonTrait, g1, ::LineStringTrait, g2)::Bool = line_crosses_poly(g2, g1)

function multipoint_crosses_line(geom1, geom2)
    int_point = false
    ext_point = false
    i = 1
    np2 = GI.npoint(geom2)

    while i < GI.npoint(geom1) && !int_point && !ext_point
        for j in 1:GI.npoint(geom2) - 1
            exclude_boundary = (j === 1 || j === np2 - 2) ? :none : :both
            if point_on_segment(GI.getpoint(geom1, i), (GI.getpoint(geom2, j), GI.getpoint(geom2, j + 1)); exclude_boundary)
                int_point = true
            else
                ext_point = true
            end
        end
        i += 1
    end

    return int_point && ext_point
end

function line_crosses_line(line1, line2)
    np2 = GI.npoint(line2)
    if GeometryOps.intersects(line1, line2)
        for i in 1:GI.npoint(line1) - 1
            for j in 1:GI.npoint(line2) - 1
                exclude_boundary = (j === 1 || j === np2 - 2) ? :none : :both
                pa = GI.getpoint(line1, i)
                pb = GI.getpoint(line1, i + 1)
                p = GI.getpoint(line2, j)
                point_on_segment(p, (pa, pb); exclude_boundary) && return true
            end
        end
    end
    return false
end

function line_crosses_poly(line, poly)
    for l in flatten(AbstractCurveTrait, poly)
        intersects(line, l) && return true
    end
    return false
end

function multipoint_crosses_poly(mp, poly)
    int_point = false
    ext_point = false

    for p in GI.getpoint(mp)
        if point_in_polygon(p, poly)
            int_point = true
        else
            ext_point = true
        end
        int_point && ext_point && return true
    end
    return false
end
