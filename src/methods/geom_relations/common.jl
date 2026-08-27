# Identical boilerplate methods for geom relations live here

#= The manifold-threaded family: every forwarding method carries the manifold
through to the underlying call. `crosses` and `overlaps` are generated below
instead, without a manifold, because they still run on their own planar helpers
rather than the shared processors. =#
for f in (:coveredby, :disjoint, :touches, :within)
    _f = Symbol(:_, f)

    @eval begin
        # Features
        $_f(m::Manifold, ::GI.FeatureTrait, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(m, GI.geometry(g1), g2; kw...)
        $_f(m::Manifold, ::GI.AbstractGeometryTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(m, g1, GI.geometry(g2); kw...)
        $_f(m::Manifold, ::GI.FeatureTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(m, GI.geometry(g1), GI.geometry(g2); kw...)

        # Extent forwarding
        $_f(m::Manifold, t1::GI.FeatureTrait, f1, ::GI.RectangleTrait, e::Extents.Extent; kw...) =
            $_f(m, t1, f1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::GI.RectangleTrait, e1::Extents.Extent, t2::GI.FeatureTrait, f2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, f2; kw...)
        $_f(m::Manifold, t1::GI.AbstractGeometryTrait, g1, ::GI.RectangleTrait, e::Extents.Extent; kw...) =
            $_f(m, t1, g1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::GI.RectangleTrait, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, g2; kw...)
        #= Two bare extents carry no manifold information — they are coordinate
        boxes, not geometries — so both manifolds answer with the box predicate. =#
        $_f(m::Manifold, ::GI.RectangleTrait, e1::Extents.Extent, ::GI.RectangleTrait, e2::Extents.Extent; kw...) =
            Extents.$f(e1, e2)

        # Backwards compatibility for when Extent traits were Nothing
        $_f(m::Manifold, t1::GI.FeatureTrait, f1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(m, t1, f1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::Nothing, e1::Extents.Extent, t2::GI.FeatureTrait, f2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, f2; kw...)
        $_f(m::Manifold, t1::GI.AbstractGeometryTrait, g1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(m, t1, g1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::Nothing, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, g2; kw...)
        $_f(m::Manifold, ::Nothing, e1::Extents.Extent, ::Nothing, e2::Extents.Extent; kw...) =
            Extents.$f(e1, e2)

        # Table rows ? or error
        $_f(m::Manifold, ::Nothing, g1, ::GI.FeatureTrait, f2; kw...) = $f(m, _geometry_or_error(g1; kw...), f2)
        $_f(m::Manifold, ::GI.FeatureTrait, f1, ::Nothing, g2; kw...) = $f(m, f1, _geometry_or_error(g2; kw...))
        $_f(m::Manifold, ::Nothing, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(m, _geometry_or_error(g1; kw...), g2)
        $_f(m::Manifold, ::GI.AbstractGeometryTrait, g1, ::Nothing, g2; kw...) = $f(m, g1, _geometry_or_error(g2; kw...))
        $_f(m::Manifold, ::Nothing, g1, ::Nothing, g2; kw...) =
            $f(m, _geometry_or_error(g1; kw...), _geometry_or_error(g2; kw...))
    end
end

for f in (:crosses, :overlaps)
    _f = Symbol(:_, f)

    @eval begin
        # Features
        $_f(::GI.FeatureTrait, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(GI.geometry(g1), g2; kw...)
        $_f(::GI.AbstractGeometryTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(g1, GI.geometry(g2); kw...)
        $_f(::GI.FeatureTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(GI.geometry(g1), GI.geometry(g2); kw...)

        # Extent forwarding
        $_f(t1::GI.FeatureTrait, f1, ::GI.RectangleTrait, e::Extents.Extent; kw...) =
            $_f(t1, f1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(::GI.RectangleTrait, e1::Extents.Extent, t2::GI.FeatureTrait, f2; kw...) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, f2; kw...)
        $_f(t1::GI.AbstractGeometryTrait, g1, ::GI.RectangleTrait, e::Extents.Extent; kw...) =
            $_f(t1, g1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(::GI.RectangleTrait, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2; kw...) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, g2; kw...)
        $_f(::GI.RectangleTrait, e1::Extents.Extent, ::GI.RectangleTrait, e2::Extents.Extent; kw...) =
            Extents.$f(e1, e2)

        # Backwards compatibility for when Extent traits were Nothing
        $_f(t1::GI.FeatureTrait, f1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(t1, f1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(::Nothing, e1::Extents.Extent, t2::GI.FeatureTrait, f2; kw...) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, f2; kw...)
        $_f(t1::GI.AbstractGeometryTrait, g1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(t1, g1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(::Nothing, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2; kw...) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, g2; kw...)
        $_f(::Nothing, e1::Extents.Extent, ::Nothing, e2::Extents.Extent; kw...) =
            Extents.$f(e1, e2)

        # Table rows ? or error
        $_f(::Nothing, g1, ::GI.FeatureTrait, f2; kw...) = $f(_geometry_or_error(g1; kw...), f2)
        $_f(::GI.FeatureTrait, f1, ::Nothing, g2; kw...) = $f(f1, _geometry_or_error(g2; kw...))
        $_f(::Nothing, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(_geometry_or_error(g1; kw...), g2)
        $_f(::GI.AbstractGeometryTrait, g1, ::Nothing, g2; kw...) = $f(g1, _geometry_or_error(g2; kw...))
        $_f(::Nothing, g1, ::Nothing, g2; kw...) =
            $f(_geometry_or_error(g1; kw...), _geometry_or_error(g2; kw...))
    end
end

@noinline _throw_no_manifold_predicate(f, m) = throw(ArgumentError(
    "`$(f)` has no $(typeof(m).name.name) implementation: it is still on the " *
    "legacy planar path rather than the shared geometry-relation processors. " *
    "Use `relate_predicate(RelateNG($m), pred_$(f)(), a, b)` instead."))
