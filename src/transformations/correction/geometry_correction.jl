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

fix(geometry) = fix(GI.trait(geometry), geometry)

function fix(trait::Trait, geometry; corrections = GeometryCorrection[ClosedRing(),], kwargs...) where Trait <: GI.AbstractGeometryTrait
    traits = TraitTarget.(application_level.(corrections))
    final_geometry = geometry
    for correction in corrections
        if trait in TraitTarget(application_level(corrections))
            final_geometry = apply()
        end
    end
    return final_geometry
end

# The API application_level exists, so from that we need to derive, given a geometry with a trait, which corrections we can apply to it.
# We can do this by running through all subtypes (recursively) of 

function _get_subtypes!(vec::Vector{DataType}, type::Type)
    if isabstracttype(type)
        for subtype in subtypes(type)
            _get_subtypes!(vec, subtype)
        end
    else # is a concrete type
        push!(vec, type)
    end
end

function _get_subtypes(type::Type)
    v = Vector{DataType}()
    _get_subtypes!(v, type)
    return v
end


# ## Available corrections

#=
```@autodocs; canonical = false
Modules = [GeometryOps]
Filter = t -> typeof(t) === DataType && t <: GeometryOps.GeometryCorrection
```
=#

#=

Old code:

This code was meant to batch corrections across multiple geometries. However, this
fails when it encounters things like TraitTargets across multiple levels.


=#