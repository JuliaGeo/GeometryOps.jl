# Identical boilerplate methods for geom relations live here

for f in (:coveredby, :crosses, :disjoint, :overlaps, :touches, :within)
    _f = Symbol(:_, f)

    @eval begin
        # Features (4-arg form without manifold)
        $_f(::GI.FeatureTrait, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(GI.geometry(g1), g2; kw...)
        $_f(::GI.AbstractGeometryTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(g1, GI.geometry(g2); kw...)
        $_f(::GI.FeatureTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(GI.geometry(g1), GI.geometry(g2); kw...)

        # Features (5-arg form with manifold)
        $_f(m::Manifold, ::GI.FeatureTrait, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(m, GI.geometry(g1), g2; kw...)
        $_f(m::Manifold, ::GI.AbstractGeometryTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(m, g1, GI.geometry(g2); kw...)
        $_f(m::Manifold, ::GI.FeatureTrait, g1, ::GI.FeatureTrait, g2; kw...) = $f(m, GI.geometry(g1), GI.geometry(g2); kw...)

        # Extent forwarding (4-arg form without manifold)
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

        # Extent forwarding (5-arg form with manifold)
        $_f(m::Manifold, t1::GI.FeatureTrait, f1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(m, t1, f1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::Nothing, e1::Extents.Extent, t2::GI.FeatureTrait, f2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, f2; kw...)
        $_f(m::Manifold, t1::GI.AbstractGeometryTrait, g1, ::Nothing, e::Extents.Extent; kw...) =
            $_f(m, t1, g1, GI.PolygonTrait(), extent_to_polygon(e); kw...)
        $_f(m::Manifold, ::Nothing, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2; kw...) =
            $_f(m, GI.PolygonTrait(), extent_to_polygon(e1), t2, g2; kw...)
        $_f(::Manifold, ::Nothing, e1::Extents.Extent, ::Nothing, e2::Extents.Extent; kw...) =
            Extents.$f(e1, e2)

        # Table rows ? or error (4-arg form)
        $_f(::Nothing, g1, ::GI.FeatureTrait, f2; kw...) = $f(_geometry_or_error(g1; kw...), f2)
        $_f(::GI.FeatureTrait, f1, ::Nothing, g2; kw...) = $f(f1, _geometry_or_error(g2; kw...))
        $_f(::Nothing, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(_geometry_or_error(g1; kw...), g2)
        $_f(::GI.AbstractGeometryTrait, g1, ::Nothing, g2; kw...) = $f(g1, _geometry_or_error(g2; kw...))
        $_f(::Nothing, g1, ::Nothing, g2; kw...) =
            $f(_geometry_or_error(g1; kw...), _geometry_or_error(g2; kw...))

        # Table rows ? or error (5-arg form with manifold)
        $_f(m::Manifold, ::Nothing, g1, ::GI.FeatureTrait, f2; kw...) = $f(m, _geometry_or_error(g1; kw...), f2)
        $_f(m::Manifold, ::GI.FeatureTrait, f1, ::Nothing, g2; kw...) = $f(m, f1, _geometry_or_error(g2; kw...))
        $_f(m::Manifold, ::Nothing, g1, ::GI.AbstractGeometryTrait, g2; kw...) = $f(m, _geometry_or_error(g1; kw...), g2)
        $_f(m::Manifold, ::GI.AbstractGeometryTrait, g1, ::Nothing, g2; kw...) = $f(m, g1, _geometry_or_error(g2; kw...))
        $_f(m::Manifold, ::Nothing, g1, ::Nothing, g2; kw...) =
            $f(m, _geometry_or_error(g1; kw...), _geometry_or_error(g2; kw...))
        end
end
