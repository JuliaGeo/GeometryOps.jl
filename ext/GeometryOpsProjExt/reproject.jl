import GeometryOps: GI, GeoInterface, reproject, apply, transform, _is3d, True, False, booltype, ThreadFunctors
import GeoFormatTypes
import Proj

# TODO:
# - respect `time`
# - respect measured values

struct ApplyToPoint{Z, F}
    f::F
end

ApplyToPoint{Z}(f::F) where {Z, F} = ApplyToPoint{Z, F}(f)

(t::ApplyToPoint{false})(p) = t.f(GI.x(p), GI.y(p))
(t::ApplyToPoint{true})(p) = t.f(GI.x(p), GI.y(p), GI.z(p))


function reproject(geom;
    source_crs=nothing, target_crs=nothing, transform=nothing, kw...
)
    if isnothing(transform)
        if isnothing(source_crs) 
            source_crs = if GI.trait(geom) isa Nothing && geom isa AbstractArray
                GeoInterface.crs(first(geom))
            else
                GeoInterface.crs(geom)
            end
        end

        # If its still nothing, error
        isnothing(source_crs) && throw(ArgumentError("geom has no crs attached. Pass a `source_crs` keyword"))

        # Otherwise reproject
        reproject(geom, source_crs, target_crs; kw...)
    else
        reproject(geom, transform; kw...)
    end
end

function reproject(geom, source_crs, target_crs;
    time=Inf,
    threaded=False(),
    always_xy=true,
    kw...
)
    return reproject(geom, Proj.Transformation(source_crs, target_crs; always_xy); target_crs, time, threaded, kw...)
end

function reproject(geom, transform::Proj.Transformation; target_crs = nothing, time=Inf, threaded = False(), kw...)
    if isnothing(target_crs)
        target_crs = GeoFormatTypes.ESRIWellKnownText(Proj.CRS(Proj.proj_get_target_crs(transform.pj)))
    end
    if booltype(threaded) isa True
        isnothing(transform) || throw(ArgumentError("threaded reproject doesn't accept a single Transformation"))
        tasks_per_thread = 2
        ntasks = Threads.nthreads() * tasks_per_thread
        functors = [ApplyToPoint{_is3d(geom)}(Proj.Transformation(Proj.proj_clone(Proj.proj_context_clone(), transform.pj))) for _ in 1:ntasks]
        transforms = ThreadFunctors(functors, tasks_per_thread)
        return apply(transforms, GI.PointTrait(), geom; crs=target_crs, kw...)
    else # threaded isa False
        return apply(ApplyToPoint{_is3d(geom)}(transform), GI.PointTrait(), geom; crs = target_crs, kw...)
    end    
end
