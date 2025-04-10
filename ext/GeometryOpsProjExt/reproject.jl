import GeometryOps: GI, GeoInterface, reproject, apply, transform, _is3d, True, False, ThreadFunctors
import Proj

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
struct ApplyToPoint2{F}
    f::F
end
struct ApplyToPoint3{F}
    f::F
end

(t::ApplyToPoint2)(p) = t.f(GI.x(p), GI.y(p))
(t::ApplyToPoint3)(p) = t.f(GI.x(p), GI.y(p), GI.z(p))

function reproject(geom, source_crs, target_crs;
    time=Inf,
    threaded=False(),
    always_xy=true,
    transform=nothing,
    kw...
)
    if istrue(threaded)
        isnothing(transform) || throw(ArgumentError("threaded reproject doesn't accept a single Transformation"))
        tasks_per_thread = 2
        n = Threads.nthreads() * tasks_per_thread
        if _is3d(geom)
            functors = [ApplyToPoint3(Proj.Transformation(s, t; always_xy)) for _ in 1:n]
            transforms = ThreadFunctors(functors, tasks_per_thread)
            return apply(transforms, GI.PointTrait(), geom; crs=target_crs, kw...)
        else
            functors = [ApplyToPoint2(Proj.Transformation(s, t; always_xy)) for _ in 1:n]
            transforms = ThreadFunctors(functors, tasks_per_thread)
            return apply(transforms, GI.PointTrait(), geom; crs=target_crs, kw...)
        end
    else
        transform = if isnothing(transform) 
            s = source_crs isa Proj.CRS ? source_crs : convert(String, source_crs)
            t = target_crs isa Proj.CRS ? target_crs : convert(String, target_crs)
            Proj.Transformation(s, t; always_xy)
        else
            transform
        end
        return reproject(geom, transform; time, target_crs, kw...)
    end
end
function reproject(geom, transform::Proj.Transformation; time=Inf, target_crs=nothing, kw...)
    if _is3d(geom)
        return apply(ApplyToPoint3(transform), GI.PointTrait(), geom; crs=target_crs, kw...)
    else
        return apply(ApplyToPoint2(transform), GI.PointTrait(), geom; crs=target_crs, kw...)
    end
end
