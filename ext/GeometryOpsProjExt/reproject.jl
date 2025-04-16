import GeometryOps: GI, GeoInterface, reproject, apply, transform, _is3d, istrue,
    True, False, TaskFunctors, WithXY, WithXYZ
import GeoFormatTypes
import Proj

# TODO:
# - respect `time`
# - respect measured values

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
function reproject(geom, source_crs, target_crs; always_xy=true, kw...)
    transform = Proj.Transformation(convert(String, source_crs), convert(String, target_crs); always_xy)
    return reproject(geom, transform; target_crs, kw...)
end
function reproject(
    geom, source_crs::CRSType, target_crs::CRSType; always_xy=true, kw...
) where CRSType <: Union{GeoFormatTypes.GeoFormat, Proj.CRS}
    transform = Proj.Transformation(source_crs, target_crs; always_xy)
    return reproject(geom, transform; target_crs, kw...)
end
function reproject(geom, transform::Proj.Transformation; 
    context=C_NULL, 
    target_crs=nothing, 
    time=Inf, 
    threaded=False(), 
    kw...
)
    if isnothing(target_crs)
        target_crs = GeoFormatTypes.ESRIWellKnownText(Proj.CRS(Proj.proj_get_target_crs(transform.pj)))
    end
    kw1 = (; crs=target_crs, threaded, kw...)
    if istrue(threaded)
        tasks_per_thread = 2
        ntasks = Threads.nthreads() * tasks_per_thread
        # Clone the transformation once for each task.  
        # Currently, these transformations live in the same context, but we will soon
        # assign them to per-task contexts.
        proj_transforms = [Proj.Transformation(Proj.proj_clone(transform.pj)) for _ in 1:ntasks]

        # Construct one context per planned task
        contexts = [Proj.proj_context_clone(context) for _ in 1:ntasks]
        # Assign the context to the transformation.  We use `foreach` here 
        # to avoid generating output where we don't have to.
        foreach(Proj.proj_assign_context, getproperty.(proj_transforms, :pj), contexts)

        results = if _is3d(geom)
            functors = TaskFunctors(WithXYZ.(proj_transforms))
            apply(functors, GI.PointTrait(), geom; kw1...)
        else
            functors = TaskFunctors(WithXY.(proj_transforms))
            apply(functors, GI.PointTrait(), geom; kw1...)
        end
        # Destroy the temporary threading contexts that we created
        foreach(Proj.proj_context_destroy, contexts)
        # Return the results
        return results
    else
        if _is3d(geom)
            return apply(WithXYZ(transform), GI.PointTrait(), geom; kw1...)
        else
            return apply(WithXY(transform), GI.PointTrait(), geom; kw1...)
        end
    end    
end