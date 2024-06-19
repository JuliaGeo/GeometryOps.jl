#=
# Buffer

Buffering a geometry means computing the region `distance` away from it, and returning
that region as the new geometry.

As of now, we only support `GEOS` as the backend, meaning that LibGEOS must be loaded.
=#

function buffer(geometry, distance; kwargs...)
    buffered = buffer(GEOS(; kwargs...), geometry, distance)
    return tuples(buffered)
end

# Below is an error handler similar to the others we have for e.g. segmentize,
# which checks if there is a method error for the geos backend.


# Add an error hint for GeodesicSegments if Proj is not loaded!
function _buffer_error_hinter(io, exc, argtypes, kwargs)
    if isnothing(Base.get_extension(GeometryOps, :GeometryOpsLibGEOSExt)) && exc.f == buffer && first(argtypes) == GEOS
        print(io, "\n\nThe `buffer` method requires the LibGEOS.jl package to be explicitly loaded.\n")
        print(io, "You can do this by simply typing ")
        printstyled(io, "using LibGEOS"; color = :cyan, bold = true)
        println(io, " in your REPL, \nor otherwise loading LibGEOS.jl via using or import.")
    end
end
