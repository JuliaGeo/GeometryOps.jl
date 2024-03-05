# # Geometry Corrections

export fix

# This file simply defines the `GeometryCorrection` abstract type, and the interface that any `GeometryCorrection` must implement.
#=

A geometry correction is a transformation that is applied to a geometry to correct it in some way. 

For example, a [`ClosedRing`](@ref) correction might be applied to a `Polygon` to ensure that its exterior ring is closed.

## Interface

All `GeometryCorrection`s are callable structs which, when called, apply the correction to the given geometry, and return either a copy or the original geometry (if nothing needed to be corrected).

See below for the full interface specification.

```@docs; canonical = false
GeometryOps.GeometryCorrection
```

Any geometry correction must implement the interface as given above. 



=#

"""
    abstract type GeometryCorrection

This abstract type represents a geometry correction.

## Interface

Any `GeometryCorrection` must implement two functions:
    * `application_level(::GeometryCorrection)::AbstractGeometryTrait`: This function should return the `GeoInterface` trait that the correction is intended to be applied to, like `PointTrait` or `LineStringTrait` or `PolygonTrait`.
    * `(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)`: This function should apply the correction to the given geometry, and return a new geometry.
"""
abstract type GeometryCorrection end

application_level(gc::GeometryCorrection) = error("Not implemented yet for $(gc)")

(gc::GeometryCorrection)(geometry) = gc(GI.trait(geometry), geometry)

(gc::GeometryCorrection)(trait::GI.AbstractGeometryTrait, geometry) = error("Not implemented yet for $(gc) and $(trait).")

function fix(geometry; corrections = GeometryCorrection[ClosedRing(),], kwargs...)
    traits = application_level.(corrections)
    final_geometry = geometry
    for Trait in (GI.PointTrait, GI.MultiPointTrait, GI.LineStringTrait, GI.LinearRingTrait, GI.MultiLineStringTrait, GI.PolygonTrait, GI.MultiPolygonTrait)
        available_corrections = findall(x -> x == Trait, traits)
        isempty(available_corrections) && continue
        println("Correcting for $(Trait)")
        net_function = reduce(∘, corrections[available_corrections])
        final_geometry = apply(net_function, Trait, final_geometry; kwargs...)
    end
    return final_geometry
end

# ## Available corrections

#=
```@autodocs; canonical = false
Modules = [GeometryOps]
Filter = t -> typeof(t) === DataType && t <: GeometryOps.GeometryCorrection
```
=#