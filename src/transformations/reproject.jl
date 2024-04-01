# # Geometry reprojection

export reproject

# This file is pretty simple - it simply reprojects a geometry pointwise from one CRS
# to another. It uses the `Proj` package for the transformation, but this could be 
# moved to an extension if needed.

# Note that the actual implementation is in the `GeometryOpsProjExt` extension module.

# This works using the [`apply`](@ref) functionality.

"""
    reproject(geometry; source_crs, target_crs, transform, always_xy, time)
    reproject(geometry, source_crs, target_crs; always_xy, time)
    reproject(geometry, transform; always_xy, time)

Reproject any GeoInterface.jl compatible `geometry` from `source_crs` to `target_crs`.

The returned object will be constructed from `GeoInterface.WrapperGeometry`
geometries, wrapping views of a `Vector{Proj.Point{D}}`, where `D` is the dimension.

!!! tip
    The `Proj.jl` package must be loaded for this method to work, 
    since it is implemented in a package extension.

## Arguments

- `geometry`: Any GeoInterface.jl compatible geometries.
- `source_crs`: the source coordinate referece system, as a GeoFormatTypes.jl object or a string.
- `target_crs`: the target coordinate referece system, as a GeoFormatTypes.jl object or a string.

If these a passed as keywords, `transform` will take priority.
Without it `target_crs` is always needed, and `source_crs` is
needed if it is not retreivable from the geometry with `GeoInterface.crs(geometry)`.

## Keywords

- `always_xy`: force x, y coordinate order, `true` by default.
    `false` will expect and return points in the crs coordinate order.
- `time`: the time for the coordinates. `Inf` by default.
$APPLY_KEYWORDS
"""
function reproject end

# ## Method error handling
# We also inject a method error handler, which 
# prints a suggestion if the Proj extension is not loaded.

function _reproject_error_hinter(io, exc, argtypes, kwargs)
    if isnothing(Base.get_extension(GeometryOps, :GeometryOpsProjExt)) && exc.f == reproject
        print(io, "\n\nThe `reproject` method requires the Proj.jl package to be explicitly loaded.\n")
        print(io, "You can do this by simply typing ")
        printstyled(io, "using Proj"; color = :cyan, bold = true)
        println(io, " in your REPL, \nor otherwise loading Proj.jl via using or import.")
    else # this is a more general error
        nothing
    end
end