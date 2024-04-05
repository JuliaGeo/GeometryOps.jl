#=
# Ring orientation

This file contains a test for linear ring orientation in polygons.  Specifically, the exterior should be clockwise and the interior should be anticlockwise.

=#


"""
    RingOrientation() <: GeometryCorrection

This correction ensures that a polygon's exterior ring is oriented clockwise, and its interior rings are all oriented 

It can be called on any geometry correction as usual.

See also [`GeometryCorrection`](@ref).
"""
struct RingOrientation <: GeometryCorrection end

application_level(::RingOrientation) = GI.PolygonTrait()

function _check_and_reverse_hole(hole)
    if isclockwise(hole)
        return GI.LinearRing(reverse(GI.getpoint(hole)))
    else
        return GI.LinearRing((GI.getpoint(hole)))
    end
end

function (::RingOrientation)(::GI.PolygonTrait, polygon)
    # First, test the exterior
    exterior = GI.getexterior(polygon)
    if !isclockwise(exterior) # exterior must be clockwise
      exterior = GI.LinearRing(reverse(GI.getpoint(exterior)))
    end

    return GI.Polygon(SVector((exterior, (_check_and_reverse_hole(hole) for hole in GI.gethole(polygon))...)))
end
