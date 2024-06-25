# # PolygonContents

export PolygonContents

#=
Polygons should only contain linear rings.  This fix checks
whether the contents of the polygon are linear rings or linestrings,
and converts linestrings to linear rings.

It does **NOT** check whether the linear rings are valid - it only checks
the types of the polygon's constituent geometries.  You can use the [`ClosedRing`](@ref)
geometry fix to check for validity after applying this fix.
=#

struct PolygonContents <: GeometryCorrection end

application_level(::PolygonContents) = GI.PolygonTrait

function (::PolygonContents)(::GI.PolygonTrait, polygon)
    exterior = GI.getexterior(polygon)
    fixed_exterior = _ls2lr(exterior)
    holes = GI.gethole(polygon)
    if isempty(holes)
        return GI.Polygon([fixed_exterior])
    end
    fixed_holes = _ls2lr.(holes)
    return GI.Polygon([fixed_exterior, fixed_holes...])
end

_ls2lr(x) = _ls2lr(GI.geomtrait(x), x)

_ls2lr(::GI.LineStringTrait, x) = GI.LinearRing(GI.getpoint(x))
_ls2lr(::GI.LinearRingTrait, x) = x