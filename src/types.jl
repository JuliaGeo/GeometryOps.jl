#=
# Types

This file defines some fundamental types used in GeometryOps.

!!! warning
    Unlike in other Julia packages, only some types are defined in this file, not all. 
    This is because we define types in the files where they are used, to make it easier to understand the code.

=#
export GEOS
#=

## `GEOS`

`GEOS` is a struct which instructs the method it's passed to as an algorithm
to use the appropriate GEOS function via `LibGEOS.jl` for the operation.

It's generally a lot slower than the native Julia implementations, but it's
useful for two reasons:
1. Functionality which doesn't exist in GeometryOps can be accessed through the GeometryOps API, but use GEOS in the backend until someone implements a native Julia version.
2. It's a good way to test the correctness of the native implementations.

=#

"""
    GEOS(; params...)

A struct which instructs the method it's passed to as an algorithm
to use the appropriate GEOS function via `LibGEOS.jl` for the operation.

Dispatch is generally carried out using the names of the keyword arguments.
For example, `segmentize` will only accept a `GEOS` struct with only a 
`max_distance` keyword, and no other.

It's generally a lot slower than the native Julia implementations, since
it must convert to the LibGEOS implementation and back - so be warned!
"""
struct GEOS
    params::NamedTuple
end

function GEOS(; params...)
    nt = NamedTuple(params)
    return GEOS(nt)
end
# These are definitions for convenience, so we don't have to type out 
# `alg.params` every time.
Base.get(alg::GEOS, key, value) = Base.get(alg.params, key, value)
Base.get(f::Function, alg::GEOS, key) = Base.get(f, alg.params, key)

"""
    enforce(alg::GO.GEOS, kw::Symbol, f)

Enforce the presence of a keyword argument in a `GEOS` algorithm, and return `alg.params[kw]`.

Throws an error if the key is not present, and mentions `f` in the error message (since there isn't 
a good way to get the name of the function that called this method).
"""
function enforce(alg::GEOS, kw::Symbol, f)
    if haskey(alg.params, kw)
        return alg.params[kw]
    else
        error("$(f) requires a `$(kw)` keyword argument to the `GEOS` algorithm, which was not provided.")
    end
end
