# Identical boilerplate methods for geom relations live here

for f in (:coveredby, :crosses, :disjoint, :overlaps, :touches, :within)
    _f = Symbol(:_, f)

    @eval begin
        # Features
        $_f(::GI.FeatureTrait, g1, ::GI.AbstractGeometryTrait, g2) = $f(geometry(g1), g2)
        $_f(::GI.AbstractGeometryTrait, g1, ::GI.FeatureTrait, g2) = $f(g1, geometry(g2))
        $_f(::GI.FeatureTrait, g1, ::GI.FeatureTrait, g2) = $f(geometry(g1), geometry(g2))

        # Extent forwarding
        $_f(t1::GI.FeatureTrait, f1, ::Nothing, e::Extents.Extent) =
            $_f(t1, f1, GI.PolygonTrait(), extent_to_polygon(e))
        $_f(::Nothing, e1::Extents.Extent, t2::GI.FeatureTrait, f2) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, f2)
        $_f(t1::GI.AbstractGeometryTrait, g1, ::Nothing, e::Extents.Extent) =
            $_f(t1, g1, GI.PolygonTrait(), extent_to_polygon(e))
        $_f(::Nothing, e1::Extents.Extent, t2::GI.AbstractGeometryTrait, g2) =
            $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, g2)
        $_f(::Nothing, e1::Extents.Extent, ::Nothing, e2::Extents.Extent) =
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
