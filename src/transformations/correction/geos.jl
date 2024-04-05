#=
# LibGEOS correction

This file defines exportable skeleton code for a "GEOSCorrection" which 
simply invokes `LibGEOS.makeValid` on the geometry.

The glue code is defined in the extension on LibGEOS.

=#


"""
    GEOSCorrection() <: GeometryCorrection

This correction runs `LibGEOS.makeValid` on the highest level of each geometry available.

See also [`GeometryCorrection`](@ref).
"""
struct GEOSCorrection <: GeometryCorrection end

application_level(::GEOSCorrection) = TraitTarget(GI.MultiPolygonTrait(), GI.PolygonTrait(), GI.LinearRingTrait(), GI.MultiLineStringTrait(), GI.LineStringTrait(), GI.MultiPointTrait(), GI.PointTrait())

