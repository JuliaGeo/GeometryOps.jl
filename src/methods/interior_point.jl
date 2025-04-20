struct YSpanInteriorPoint{M <: Manifold} <: Algorithm{M}
    manifold::M
    # accelerator::IntersectionAccelerator # from clipping_processor
    # selection_criterion::Symbol
end

YSpanInteriorPoint() = YSpanInteriorPoint(Planar(), AutoAccelerator())

manifold(a::YSpanInteriorPoint) = a.manifold
best_manifold(a::YSpanInteriorPoint, input) = Planar()
rebuild(a::YSpanInteriorPoint, manifold::Manifold) = YSpanInteriorPoint(manifold, a.accelerator)
rebuild(a::YSpanInteriorPoint, accelerator::IntersectionAccelerator) = YSpanInteriorPoint(a.manifold, accelerator)

interior_point(geom, kw...) = interior_point(YSpanInteriorPoint(), geom; kw...)
function interior_point(alg::Algorithm, geom; kw...)
    return apply(GeometryOpsCore.ApplyKnowingTrait(Base.Fix1(interior_point, alg)), TraitTarget{GI.AbstractGeometryTrait}(), geom; kw...)
end

interior_point(alg::YSpanInteriorPoint, trait::GI.AbstractGeometryTrait, geom) = error("Not implemented yet for trait $trait and algorithm $alg - file an issue on GitHub if you need this!")

interior_point(alg::YSpanInteriorPoint, ::GI.PointTrait, geom) = geom
interior_point(alg::YSpanInteriorPoint, ::GI.MultiPointTrait, geom) = rand(GI.getgeom(geom))


function __extent_intersects_y(extent, y)
    return extent.Y[1] <= y <= extent.Y[2]
end

#=
For polygons, this algorithm actually performs two passes over the edges.
First, it finds the two y coordinates that are closest to the center of the polygon's 
y-extent, then it takes the average of those to get the "representative" Y coordinate.

Once that is done, it runs a scanline algorithm to find all the edges that intersect
the line at the Y coordinate it found.  All these edges are then processed and the 
point that is returned is the midpoint of the longest line segment of the scanline
that is within the polygon.
=#
function interior_point(alg::YSpanInteriorPoint{Planar}, ::Union{GI.PolygonTrait, GI.LinearRingTrait}, geom)
    geom_extent = GI.extent(geom)
    best_y = _find_best_y_point(geom)

    intersecting_edges, _intersecting_idxs = to_edgelist(Extents.Extent(X = geom_extent.X, Y = (best_y, best_y)), geom)
    intersection_xs = map(intersecting_edges) do edge
        _interpolate_x_to_line(GI.getgeom(edge, 1), GI.getgeom(edge, 2), best_y)
    end

    sort!(intersection_xs)

    max_width = -Inf
    max_idx = 0

    for i in 1:2:length(intersection_xs)
        width = abs(intersection_xs[i+1] - intersection_xs[i])
        if width > max_width
            max_width = width
            max_idx = i
        end
    end

    if max_idx == 0
        error("Weird thing happened in `intersection_point`, max_width is 0")
    end

    return ((intersection_xs[max_idx] + intersection_xs[max_idx+1]) / 2, best_y)
end

function _find_best_y_point(geom)
    geom_extent = GI.extent(geom)
    ylo, yhi = extrema(geom_extent.Y)
    ycenter = (ylo + yhi) / 2

    for point in GI.getpoint(geom)
        current_y = GI.y(point)
        if current_y <= ycenter
            if current_y > ylo
                ylo = current_y
            end
        else # current_y > ycenter
            if current_y < yhi
                yhi = current_y
            end
        end
    end

    return (ylo + yhi) / 2
    
end

function _interpolate_x_to_line(p1, p2, y)
    x0 = GI.x(p1)
    x1 = GI.x(p2)

    if x0 == x1
        return x0
    end
    
    # Assert: segDX is non-zero, due to previous equality test
    segDX = x1 - x0
    segDY = GI.y(p1) - GI.y(p2)
    m = segDY / segDX
    x = x0 + ((y - GI.y(p1)) / m)
    return x
end