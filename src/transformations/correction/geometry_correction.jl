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
  - `application_level(::GeometryCorrection)::TraitTarget`: This function should 
  return the `GeoInterface` trait that the correction is intended to be applied to, 
  like `PointTrait` or `LineStringTrait` or `PolygonTrait`.  It can also return a 
  union of traits via `TraitTarget`, but that behaviour is a bit tricky...
  - `(::GeometryCorrection)(::AbstractGeometryTrait, geometry)::(some_geometry)`: 
  This function should apply the correction to the given geometry, and return a new 
  geometry.
"""
abstract type GeometryCorrection end

# Make sure that geometry corrections are treated as scalars when broadcasting.
Base.Broadcast.broadcastable(c::GeometryCorrection) = (c,)

application_level(gc::GeometryCorrection) = error("Not implemented yet for $(gc)")

(gc::GeometryCorrection)(geometry) = gc(GI.trait(geometry), geometry)

(gc::GeometryCorrection)(trait::GI.AbstractGeometryTrait, geometry) = error("Not implemented yet for $(gc) and $(trait).")

function fix(geometry; corrections = GeometryCorrection[ClosedRing()], kwargs...)
    final_geoms = geometry
    # Iterate through the corrections and apply them to the input.
    # This allocates a _lot_, especially when reconstructing tables,
    # but it's the only fully general way to do this that I can think of.
    for correction in corrections
        final_geoms = apply(correction, application_level(correction), final_geoms; kwargs...)
    end
    #=
    # This was the old implementation
    application_levels = application_level.(corrections)
    final_geometry = geometry
    for trait in (GI.PointTrait(), GI.MultiPointTrait(), GI.LineStringTrait(), GI.LinearRingTrait(), GI.MultiLineStringTrait(), GI.PolygonTrait(), GI.MultiPolygonTrait())
        available_corrections = findall(x -> trait in x, application_levels)
        isempty(available_corrections) && continue
        @debug "Correcting for $(trait), with corrections: " available_corrections
        net_function = reduce(∘, corrections[available_corrections])
        # TODO: this allocates too much, because it keeps reconstructing higher level geoms.
        # We might want some way to embed the fixes in reconstruct/rebuild, which would imply a modified apply pipeline...
        final_geometry = apply(net_function, trait, final_geometry; kwargs...)
    end
    return final_geometry
    =#
    return final_geoms
end

# ## Available corrections

#=
```@autodocs; canonical = false
Modules = [GeometryOps]
Filter = t -> typeof(t) === DataType && t <: GeometryOps.GeometryCorrection
```
=#