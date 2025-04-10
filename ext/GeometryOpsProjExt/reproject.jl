import GeometryOps: GI, GeoInterface, reproject, apply, transform, _is3d, True, False, booltype, TaskFunctors
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
    return reproject(geom, Proj.Transformation(convert(String, source_crs), convert(String, target_crs); always_xy); target_crs, time, threaded, kw...)
end

function reproject(geom, source_crs::CRSType, target_crs::CRSType;
    time=Inf,
    threaded=False(),
    always_xy=true,
    kw...
) where CRSType <: Union{GeoFormatTypes.GeoFormat, Proj.CRS}
    return reproject(geom, Proj.Transformation(source_crs, target_crs; always_xy); target_crs, time, threaded, kw...)
end

function reproject(geom, transform::Proj.Transformation; context = C_NULL, target_crs = nothing, time=Inf, threaded = False(), kw...)
    if isnothing(target_crs)
        target_crs = GeoFormatTypes.ESRIWellKnownText(Proj.CRS(Proj.proj_get_target_crs(transform.pj)))
    end
    if booltype(threaded) isa True
        tasks_per_thread = 2
        ntasks = Threads.nthreads() * tasks_per_thread
        # Construct one context per planned task
        contexts = [Proj.proj_context_clone(context) for _ in 1:ntasks]
        # Clone the transformation for each context
        proj_transforms = [Proj.Transformation(Proj.proj_clone(transform.pj)) for context in contexts]
        # Assign the context to the transformation
        Proj.proj_assign_context.(getproperty.(proj_transforms, :pj), contexts)

        appliers = if _is3d(geom)
            ApplyToPoint{true}.(proj_transforms)
        else
            ApplyToPoint{false}.(proj_transforms)
        end

        functors = TaskFunctors(appliers)
        results = apply(functors, GI.PointTrait(), geom; crs=target_crs, threaded, kw...)
        # Destroy the temporary threading contexts that we created
        Proj.proj_destroy.(contexts)
        # Return the results
        return results
    else # threaded isa False
        applier = if _is3d(geom)
            ApplyToPoint{true}(transform)
        else
            ApplyToPoint{false}(transform)
        end
        return apply(applier, GI.PointTrait(), geom; threaded, crs = target_crs, kw...)
    end    
end
