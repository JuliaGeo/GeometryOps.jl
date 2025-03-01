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
# TODO: Add working example
```
"""
crosses(g1, g2)::Bool = crosses(trait(g1), g1, trait(g2), g2)::Bool
crosses(g1) = Base.Fix1(crosses, g1)

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
            if _point_on_segment(GI.getpoint(geom1, i), (GI.getpoint(geom2, j), GI.getpoint(geom2, j + 1)); exclude_boundary)
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
                _point_on_segment(p, (pa, pb); exclude_boundary) && return true
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
        if _point_polygon_process(
            p, poly;
            in_allow = true, on_allow = true, out_allow = false, exact = False()
        )
            int_point = true
        else
            ext_point = true
        end
        int_point && ext_point && return true
    end
    return false
end

#= TODO: Once crosses is swapped over to use the geom relations workflow, can
delete these helpers. =#

function _point_on_segment(point, (start, stop); exclude_boundary::Symbol=:none)::Bool
    x, y = GI.x(point), GI.y(point)
    x1, y1 = GI.x(start), GI.y(start)
    x2, y2 = GI.x(stop), GI.y(stop)

    dxc = x - x1
    dyc = y - y1
    dx1 = x2 - x1
    dy1 = y2 - y1

    # TODO use better predicate for crossing here
    cross = dxc * dy1 - dyc * dx1
    cross != 0 && return false

    # Will constprop optimise these away?
    if exclude_boundary === :none
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x <= x2 : x2 <= x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y <= y2 : y2 <= y && y <= y1
    elseif exclude_boundary === :start
        if abs(dx1) >= abs(dy1)
             return dx1 > 0 ? x1 < x && x <= x2 : x2 <= x && x < x1
        end
        return dy1 > 0 ? y1 < y && y <= y2 : y2 <= y && y < y1
    elseif exclude_boundary === :end
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x < x2 : x2 < x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y < y2 : y2 < y && y <= y1
    elseif exclude_boundary === :both
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 < x && x < x2 : x2 < x && x < x1
        end
        return dy1 > 0 ? y1 < y && y < y2 : y2 < y && y < y1
    end
    return false
end
