# Developer documentation

This is mostly some notes about how GeometryOps is structured, geared towards developers 
and folks interested in learning about the internals of GeometryOps.

## Coding standards

Every source file is also compiled into the documentation via Literate.jl, so 
it should have a Markdown section at the top that describes what that file does.  

For most cases, each function has a single file, and so that file should:
- have a header which is the capitalized name of that function
- have a docstring that is visible when you first load the page (collapsed or opened is up to you)
- have some visual examples, rendered in the documentation, of what that function does!  
  Since this is geometry, it's very easy to plot with e.g Makie.

We also request that you define [`Algorithm`](@ref) types and use those to define the behaviour of your function.  We don't have an operations interface yet (coming soon!) but that should be done as well!

## Geometry representation

In Julia there is no one fixed geometry representation; instead we can use any memory layout in a standardized way via [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl).  This means
iterating over `enumerate(GI.getpoint(geom))` rather than `1:length(geom)`, and similar. 

However, sometimes you want a fast, in-Julia geometry with known layout.  For example, `getpoint` on 
GEOS or GDAL geoms is quite slow because you have to make a `ccall` with a pointer.  For cases like these
you can simply use `GO.tuples` to convert geometries to [GeoInterface wrapper geometries](https://juliageo.org/GeoInterface.jl/dev/guides/defaults/#Wrapper-types) that wrap (usually) 2-tuple points in vectors.

When returning geometries you can generally return any GeoInterface geometry but should prefer such
GeoInterface tuple geometries as mentioned above.

