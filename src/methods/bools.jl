# These are all adapted from Turf.jl
# The may not necessarily be what want in the end but work for now

"""
    isclockwise(line::Union{LineString, Vector{Position}})::Bool

Take a ring and return true or false whether or not the ring is clockwise or counter-clockwise.

# Examples
```jldoctest
julia> line = LineString([[0, 0], [1, 1], [1, 0], [0, 0]])
LineString(Array{Float64,1}[[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]])

julia> clockwise(line)
true
```
"""
isclockwise(geom)::Bool isclockwise(GI.trait(line), geom)
function isclockwise(::AbstractCurveTrait, line)::Bool
    sum = 0.0
    prev = GI.getpoint(line, 1)
    for p in GI.getpoint(line)
        # sum will be zero for the first point as x is subtracted from itself
        sum += (GI.x(p) - GI.x(prev)) * (GI.y(p) + GI.y(prev))
        prev = p
    end

    return sum > 0.0
end

"""
    isconcave(poly::Polygon)::Bool

Take a polygon and return true or false as to whether it is concave or not.

# Examples
```jldoctest
julia> poly = Polygon([[[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]]])
Polygon(Array{Array{Float64,1},1}[[[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]]])

julia> concave(poly)
false
```
"""
function isconcave(poly)::Bool
    sign = false

    exterior = GI.getexterior(poly)
    GI.npoint(exterior) <= 4 && return false

    n = GI.npoint(exterior) - 1

    for i in 1:n
        j = ((i + 1) % n) === 0 ? 1 : (i + 1) % n
        m = ((i + 2) % n) === 0 ? 1 : (i + 2) % n

        pti = GI.gepoint(exterior, i)
        ptj = GI.gepoint(exterior, j)
        ptm = GI.gepoint(exterior, m)

        dx1 = GI.x(ptm) - GI.x(ptj)
        dy1 = GI.y(ptm) - GI.y(ptj)
        dx2 = GI.x(pti) - GI.x(ptj)
        dy2 = GI.y(pti) - GI.y(ptj)

        cross = (dx1 * dy2) - (dy1 * dx2)

        if i === 0
            sign = cross > 0
        elseif sign !== (cross > 0)
            return true
        end
    end

    return false
end


function equals(geo1, geo2)
    GI.geomtrait(geo1) !== GI.geomtrait(geo2) && return false

    GI.geomtrait(geo1) isa PointTrait && return compare_points(geo1, geo2)
    GI.geomtrait(geo1) isa LineStringTrait && return compare_lines(geo1, geo2)

    error("Cant compare $(GI.trait(geo1)) and $(GI.trait(geo2)) yet")
end

function compare_points(p1, p2)
    length(p1) !== length(p2) && return false

    for i in eachindex(p1)
        round(p1[i]; digits=10) !== round(p2[i]; digits=10) && return false
    end

    return true
end

function compare_lines(p1::Vector, p2::Vector)
    # TODO: complete this
    length(p1[1]) !== length(p2[1]) && return false
end

"""
    parallel(line1::LineString, line2::LineString)::Bool

Return `true` if each segment of `line1` is parallel to the correspondent segment of `line2`

# Examples
```jldoctest
julia> line1 = LineString([[9.170356, 45.477985], [9.164434, 45.482551], [9.166644, 45.484003]])
LineString(Array{Float64,1}[[9.17036, 45.478], [9.16443, 45.4826], [9.16664, 45.484]])

julia> line2 = LineString([[9.169356, 45.477985], [9.163434, 45.482551], [9.165644, 45.484003]])
LineString(Array{Float64,1}[[9.16936, 45.478], [9.16343, 45.4826], [9.16564, 45.484]])

julia> parallel(line1, line2)
true
```
"""
function parallel(line1, line2)::Bool
    seg1 = linesegment(line1)
    seg2 = linesegment(line2)

    for i in eachindex(seg1)
        coors2 = nothing
        coors1 = seg1[i]
        coors2 = seg2i]
        _isparallel(coors1, coors2) == false && return false
    end

    return true
end

@inline function _isparallel(p1, p2)
    slope1 = bearing_to_azimuth(rhumb_bearing(GI.x(p1), GI.x(p2)))
    slope2 = bearing_to_azimuth(rhumb_bearing(GI.y(p1), GI.y(p2)))

    return slope1 === slope2
end


"""
    point_on_line(point::Point, line::LineString, ignoreEndVertices::Bool=false)::Bool

Return true if a point is on a line. Accept a optional parameter to ignore the
start and end vertices of the linestring.

# Examples
```jldoctest
julia> point = Point([1,1])
Point([1.0, 1.0])

julia> line = LineString([[0, 0], [3, 3], [4, 4]])
LineString(Array{Float64,1}[[0.0, 0.0], [3.0, 3.0], [4.0, 4.0]])

julia> point_on_line(point, line)
true
```
"""
function point_on_line(point, line; ignore_end_vertices::Bool=false)::Bool
    pCoords = point.coordinates
    lCoords = line.coordinates

    ignore = "none"
    for i in 1:length(lCoords) - 1
        if ignore_end_vertices == true
            i === 1 && (ignore = "start")
            i === length(lCoords) - 2 && (ignore = "end")
            (i === 1 && i + 1 === length(lCoords) - 1) && (ignore = "both")
        end
        point_on_segment(lCoords[i], lCoords[i + 1], pCoords, ignore) && return true
    end
    return false
end

function point_on_segment(start, stop, coords, exclude_boundary::Symbol=:none)::Bool
    x, y = GI.x(coords), GI.y(coords)
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

"""
    point_in_polygon(point::Point, polygon::Union{Polygon, MultiPolygon}, ignoreBoundary::Bool=false)::Bool

Take a Point and a Polygon and determine if the point
resides inside the polygon. The polygon can be convex or concave. The function accounts for holes.

# Examples
```jldoctest
julia> point = Point([-77, 44])
Point([-77.0, 44.0])

julia> poly = Polygon([[ [-81, 41], [-81, 47], [-72, 47], [-72, 41], [-81, 41]]])
Polygon(Array{Array{Float64,1},1}[[[-81.0, 41.0], [-81.0, 47.0], [-72.0, 47.0], [-72.0, 41.0], [-81.0, 41.0]]])

julia> point_in_polygon(point, poly)
true
```
"""
function point_in_polygon(point, polygon, ignoreBoundary::Bool=false)::Bool

    pt = point.coordinates
    poly = polygon.coordinates

    in_bbox(pt, bbox(polygon)) == false && return false

    geotype(polygon) === :Polygon && (poly = [poly])

    inside = false
    for i in eachindex(poly)
        if point_in_ring(pt, poly[i][1], ignore_boundary)
            in_hole = false
            k = 1

            while k < length(poly[i]) && !in_hole
                inRing(pt, poly[i][k], !ignore_boundary) == true && (in_hole = true)
                k += 1
            end

            !in_hole && (inside = true)
        end
    end
    return inside
end

function point_in_ring(pt, ring, ignore_boundary::Bool=false)
    inside = false

    if ring[1][1] == ring[length(ring) - 1][1] && ring[1][2] == ring[length(ring) - 1][1] 
        ring = ring[1, length(ring) - 1]
    end

    for i in 1:length(ring) - 1
        j = i + 1

        p_i = GI.getpoint(ring, i)
        p_j = GI.getpoint(ring, j)
        xi = GI.x(p_i)
        yi = GI.y(p_i)
        xj = GI.x(p_j)
        yj = GI.y(p_j)

        on_boundary = (pt[2] * (xi - xj) + yi * (xj - pt[1]) + yj * (pt[1] - xi) == 0) &&
            ((xi - pt[1]) * (xj - pt[1]) <= 0) && ((yi - pt[2]) * (yj - pt[2]) <= 0)

        on_boundary && return !ignore_boundary

        intersect = ((yi > pt[2]) !== (yj > pt[2])) && (pt[1] < (xj - xi) * (pt[2] - yi) / (yj - yi) + xi)

        intersect && (inside =  !inside)
    end

    return inside
end

function point_in_bbox(pt, bbox::Vector{Float64})
    return bbox[1] <= pt[1] &&  bbox[2] <= pt[2] &&
        bbox[3] >= pt[1] && bbox[4] >= pt[2]
end

function line_in_polygon(poly, line)
    out = false

    polybox = bbox(poly)
    linebox = bbox(line)

    !(bboxOverlap(polybox, linebox)) && return false

    coords = line.coordinates

    for i in 1:length(coords) - 1
        mid = [(coords[i][1] + coords[i + 1][1]) / 2, (coords[i][2] + coords[i + 1][2]) / 2]
        if point_in_polygon(Point(mid), poly, true)
            out = true
            break
        end
    end
    return out
end

line_on_line(tline1, line2) = line_on_line(trait(line1), line1, trait(line2), line2)
function line_on_line(t1::GI.AbstractCurveTrait, line1, t2::AbstractCurveTrait, line2)
    for p in GI.getpoint(line1)
        point_on_line(p, line2) || return false
    end
    return true
end

line_in_polygon(line, poly) = line_in_polygon(trait(line), line, trait(poly), poly)
function line_in_polygon(::LineStringTrait, line, ::PolygonTrait, poly)
    polybox = bbox(poly)
    linebox = bbox(line)

    !(bboxOverlap(polybox, linebox)) && return false

    coords = line.coordinates
    inside = false

    for i in 1:length(coords) - 1
        !(point_in_polygon(Point(coords[i]), poly)) && return false
        !inside && (inside = point_in_polygon(Point(coords[i]), poly, true))
        if !inside
            mid = [(coords[i][1] + coords[i + 1][1]) / 2, (coords[i][2] + coords[i + 1][2]) / 2]
            inside = point_in_polygon(Point(mid), poly, true)
        end
    end
    return inside
end

# TODO why were there two methods for this in Turf.jl?
function polygon_in_polygon(ft1, ft2, reverse::Bool=false)
    polybox1 = bbox(ft1)
    polybox2 = bbox(ft2)
    coords = []

    if reverse
        !(bbox_overlap(polybox2, polybox1)) && return false

        for point in GI.getpoint(ft1)
            !(point_in_polygon(point, ft2)) && return false
        end
    else
        !(bbox_overlap(polybox1, polybox2)) && return false

        for point in GI.getpoint(ft2)
            !(point_in_polygon(point, ft1)) && return false
        end
    end

    return true
end
function poly_in_poly(poly1, poly2)

    for point in GI.getpoint(poly1)
        (point_in_polygon(point, poly2)) && return true
    end

    for point in GI.getpoint(poly2)
        (point_in_polygon(point, poly1)) && return true
    end

    inter = line_intersects(polygon_to_line(poly1), polygon_to_line(poly2))
    inter != nothing && return true

    return false

end

function bbox_overlap(box1::Vector{T}, box2::Vector{T}) where {T <: Real}
    box1[1] > box2[1] && return false
    box1[3] < box2[3] && return false
    box1[2] > box2[2] && return false
    box1[4] < box2[4] && return false
    return true
end

