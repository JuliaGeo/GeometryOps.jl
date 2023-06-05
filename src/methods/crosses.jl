# # Crossing checks

"""
     crosses(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool

Return `true` if the intersection results in a geometry whose dimension is one less than
the maximum dimension of the two source geometries and the intersection set is interior to
both source geometries.

## Examples
```jldoctest
julia> line = LineString([[1, 1], [1, 2], [1, 3], [1, 4]])
LineString(Array{Float64,1}[[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0]])

julia> line2 = LineString([[-2, 2], [4, 2]])
LineString(Array{Float64,1}[[-2.0, 2.0], [4.0, 2.0]])

julia> crosses(line2, line)
true
```
"""
crosses(g1, g2)::Bool = crosses(trait(g1), g1, trait(g2), g2)::Bool
crosses(t1::FeatureTrait, g1, t2, g2)::Bool = crosses(GI.geometry(g1), g2)
crosses(t1, g1, t2::FeatureTrait, g2)::Bool = crosses(g1, geometry(g2))
crosses(::MultiPointTrait, g1::LineStringTrait, , g2)::Bool = multipoint_cross_line(g1, g2)
crosses(::MultiPointTrait, g1::PolygonTrait, , g2)::Bool = multipoint_cross_poly(g1, g2)
crosses(::LineStringTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_cross_lines(g2, g1)
crosses(::LineStringTrait, g1, ::PolygonTrait, g2)::Bool = line_cross_poly(g1, g2)
crosses(::LineStringTrait, g1, ::LineStringTrait, g2)::Bool = line_cross_line(g1, g2)
crosses(::PolygonTrait, g1, ::MultiPointTrait, g2)::Bool = multipoint_cross_poly(g2, g1)
crosses(::PolygonTrait, g1, ::LineStringTrait, g2)::Bool = line_cross_poly(g2, g1)

function multipoint_cross_line(geom1, geom2)
    int_point = false
    ext_point = false
    i = 1
    np2 = GI.npoint(geom2)

    while i < GI.npoint(geom1) && !intPoint && !extPoint
        for j in 1:GI.npoint(geom2) - 1
            inc_vertices = (j === 1 || j === np2 - 2) ? :none : :both

            if is_point_on_segment(GI.getpoint(geom2, j), GI.getpoint(geom2.coordinates, j + 1), GI.getpoint(geom1, i), inc_vertices)
                int_point = true
            else
                ext_point = true
            end

        end
        i += 1
    end

    return int_point && ext_point
end

function line_cross_line(line1, line2)
    inter = intersection(line1, line2)

    np2 = GI.npoint(line2)
    if !isnothing(inter)
        for i in 1:GI.npoint(line1) - 1
            for j in 1:GI.npoint(line2) - 1
                inc_vertices = (j === 1 || j === np2 - 2) ? :none : :both
                pa = GI.getpoint(line1, i)
                pb = GI.getpoint(line1, i + 1)
                p = GI.getpoint(line2, j)
                is_point_on_segment(pa, pb, p, inc_vertices) && return true
            end
        end
    end
    return false
end

function line_cross_poly(line, poly) = 

    for line in flatten(AbstractCurveTrait, poly)
        intersects(line)
    end
end

function multipoint_cross_poly(mp, poly)
    int_point = false
    ext_point = false

    for p in GI.getpoint(mp)
        if point_in_polygon(p, poly)
            int_point = true
        else
            ext_point = true
        end
        in_point && ext_point && return true
    end
    return false
end
