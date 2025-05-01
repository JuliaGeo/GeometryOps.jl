#=

```@meta
CollapsedDocStrings = true
```
# Applicators

Applicators are functor structs that wrap a function and a target object.  They are used to
implement the `apply` interface, and are also used to dispatch to the correct method of `_maptasks`.

The basic functor struct is `ApplyToGeom`, which just applies a function to a geometry.

```@docs; canonical=false
Applicator
ApplyToGeom
ApplyToArray
ApplyToFeatures
```
=#

"""
    abstract type Applicator{F,T}

An abstract type for applicators that apply a function to a target object.

The type parameter `F` is the type of the function to apply, and `T` is the type of the target object.

A common dispatch pattern is to dispatch on `F` which may also be e.g. a ThreadFunctor.

## Interface
All applicators must be callable by an index integer, and define the following methods:

- `rebuild(a::Applicator, f)` - swap out the function and return a new applicator.

The calling convention is `my_applicator(i::Int)`, so applicators must define this method.
"""
abstract type Applicator{F,T} end

struct ApplyToPoint{Z,M,F} <: Applicator{F,Nothing} 
    f::F
end
ApplyToPoint{Z,M}(f::F) where {Z,M,F} = ApplyToPoint{Z,M,F}(f)
ApplyToPoint{Z}(f::F) where {Z,F} = ApplyToPoint{Z,false,F}(f)
# Default function is just `tuple`
(a::Type{<:ApplyToPoint})() = a(tuple)
rebuild(a::ApplyToPoint{Z, M}, f::F) where {Z, M, F} = ApplyToPoint{Z, M, F}(f)

# Currently we ignore M by default
const WithXY = ApplyToPoint{false}
const WithXYZ = ApplyToPoint{true}
# But these could be used to require M
const WithXYM = ApplyToPoint{false,true}
const WithXYZM = ApplyToPoint{true,true}

(t::WithXY)(p) = t.f(GI.x(p), GI.y(p))
(t::WithXYZ)(p) = t.f(GI.x(p), GI.y(p), GI.z(p))
(t::WithXYZM)(p) = t.f(GI.x(p), GI.y(p), GI.z(p), GI.m(p))
(t::WithXYM)(p) = t.f(GI.x(p), GI.y(p), GI.m(p))

"""
    ApplyWithTrait(f)

ApplyWithTrait is a functor that applies a function to a trait and an object.

Specifically, the calling convention is for `f` is changed
from `f(geom)` to `f(trait, geom; kw...)`.

This is useful to keep the trait materialized through the call stack,
which can improve inferrability and performance.
"""
struct ApplyWithTrait{F} <: Applicator{F, Nothing}
    f::F
end

(a::ApplyWithTrait)(trait::GI.AbstractTrait, obj; kw...) = a.f(trait, obj; kw...)
rebuild(::ApplyWithTrait, f::F) where {F} = ApplyWithTrait{F}(f)

# ***

for T in (:ApplyToGeom, :ApplyToArray, :ApplyToFeatures, :ApplyPointsToPolygon)
    @eval begin
        struct $T{F,T,O,K} <: Applicator{F,T}
            f::F
            target::T
            obj::O
            kw::K
        end
        $T(f, target, obj; kw...) = $T(f, target, obj, kw)
        # rebuild lets us swap out the function, such as with TaskFunctors
        rebuild(a::$T, f) = $T(f, a.target, a.obj, a.kw) 
    end
end

# Functor definitions
# _maptasks may run this level threaded if `threaded==true`
# but deeper `_apply` calls will not be threaded
# For an Array there is nothing to do but map `_apply` over all values
(a::ApplyToArray)(i::Int) = _apply(a.f, a.target, a.obj[i]; a.kw..., threaded=False())
# For a FeatureCollection or Geometry we need getfeature or getgeom calls
(a::ApplyToFeatures)(i::Int) = _apply(a.f, a.target, GI.getfeature(a.obj, i); a.kw..., threaded=False())
(a::ApplyToGeom)(i::Int) = _apply(a.f, a.target, GI.getgeom(a.obj, i); a.kw..., threaded=False())
function (a::ApplyPointsToPolygon)(i::Int)
    lr = GI.getgeom(a.obj, i)
    points = map(GI.getgeom(lr)) do p
        _apply(a.f, a.target, p; a.kw..., threaded=False())
    end
    _linearring(_apply_inner(lr, points, a.kw[:crs], a.kw[:calc_extent]))
end

@doc """
    ApplyToArray(f, target, arr, kw)

Create an [`Applicator`](@ref) that applies a function to all elements of `arr`.
""" ApplyToArray

@doc """
    ApplyToGeom(f, target, geom, kw)

Create an [`Applicator`](@ref) that applies a function to all sub-geometries of `geom`.
""" ApplyToGeom

@doc """
    ApplyToFeatures(f, target, fc, kw)

Create an [`Applicator`](@ref) that applies a function to all features of `fc`.
""" ApplyToFeatures