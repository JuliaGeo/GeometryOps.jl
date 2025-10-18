# Identical boilerplate methods for geom relations live here

for f in (:coveredby, :crosses, :disjoint, :overlaps, :touches, :within)
    _f = Symbol(:_, f)

    # Features
    @eval begin
        $_f(::FeatureTrait, g1, ::Any, g2) = $f(geometry(g1), g2)
        $_f(::Any, g1, ::FeatureTrait, g2) = $f(g1, geometry(g2))
        $_f(::FeatureTrait, g1, ::FeatureTrait, g2) = $f(geometry(g1), geometry(g2))
    end

    # Table rows
    @eval begin
        $_f(::Nothing, g1, ::Any, g2) = $f(_geometry_or_error(g1), g2)
        $_f(::Any, g1, ::Nothing, g2) = $f(g1, _geometry_or_error(g2))
        $_f(::Nothing, g1, ::Nothing, g2) = 
            $f(_geometry_or_error(g1), _geometry_or_error(g2))
    end

    # Extent forwarding
    $_f(t1::GI.AbstractGeometryTrait, g1, ::Any, e::Extents.Extent) =
        $_f(t1, g1, GI.PolygonTrait(), extent_to_polygon(e))
    $_f(::Any, e1::Extents.Extent, t2::Any, g2) =
        $_f(GI.PolygonTrait(), extent_to_polygon(e1), t2, g2)
    $_f(::Any, e1::Extents.Extent, ::Any, e2::Extents.Extent) =
        Extents.$f(e1, e2)
end
