
ApplyToCoords{Z,M}(f::F) where {Z,M,F} = ApplyTo{Z,M,F}(f)
ApplyToCoords{Z}(f::F) where {Z,F} = ApplyTo{Z,false,F}(f)
# Default function is just `tuple`
(a::Type{<:ApplyToCoords}())() = a(tuple)

# Currently we ignore M by default
const ToXY = ApplyToCoords{false}
const ToXYZ = ApplyToCoords{true}
# But these could be used to require M
const ToXYM = ApplyToCoords{false,true}
const ToXYZM = ApplyToCoords{true,true}

(t::ToXY)(p) = t.f(GI.x(p), GI.y(p))
(t::ToXYZ)(p) = t.f(GI.x(p), GI.y(p), GI.z(p))
(t::ToXYZM)(p) = t.f(GI.x(p), GI.y(p), GI.m(p))
(t::ToXYM)(p) = t.f(GI.x(p), GI.y(p), GI.z(p), GI.m(p))

abstract type Applicator{F,T} end

for T in (:ApplyToGeom, :ApplyToArray, :ApplyToFeatures)
    @eval begin
        struct $T{F,T,O,K} <: Applicator{F,T}
            f::F
            target::T
            obj::O
            kw::K
        end
        $T(f, target; kw...) = $T(f, target, geom, kw)
    end
    # rebuild lets us swap out the function, such as with ThreadFunctors
    rebuild(a::Applicator, f) = $T(f, a.target, a.obj, a.kw) 
end

# Functor definitions
# _maptasks may run this level threaded if `threaded==true`
# but deeper `_apply` calls will not be threaded
# For an Array there is nothing to do but map `_apply` over all values
(a::ApplyToArray)(i::Int) = _apply(a.f, a.target, a.obj[i]; a.kw..., threaded=False())
# For a FeatureCollection or Geometry we need getfeature or getgeom calls
(a::ApplyToFeatures)(i::Int) = _apply(f, target, GI.getfeature(a.obj, i); a.kw..., threaded=False())
(a::ApplyToGeom)(i::Int) = _apply(a.f, a.target, GI.getgeom(a.obj, i); a.kw..., threaded=False())