#=
# `BoolsAsTypes`

In `apply` and `applyreduce`, we pass `threading` and `calc_extent` as types, not simple boolean values.  

This is to help compilation - with a type to hold on to, it's easier for 
the compiler to separate threaded and non-threaded code paths.

Note that if we didn't include the parent abstract type, this would have been really 
type unstable, since the compiler couldn't tell what would be returned!

We had to add the type annotation on the `_booltype(::Bool)` method for this reason as well.

TODO: should we switch to `Static.jl`?
=#

export BoolsAsTypes, True, False, booltype

"""
    abstract type BoolsAsTypes

"""
abstract type BoolsAsTypes end

"""
    struct True <: BoolsAsTypes

A struct that means `true`.
"""
struct True <: BoolsAsTypes end

"""
    struct False <: BoolsAsTypes

A struct that means `false`.
"""
struct False <: BoolsAsTypes end

# specifically for my StableTasks experiment
struct _TrueButStable <: BoolsAsTypes end

"""
    booltype(x)

Returns a [`BoolsAsTypes`](@ref) from `x`, whether it's a boolean or a BoolsAsTypes.
"""
function booltype end

@inline booltype(x::Bool)::BoolsAsTypes = x ? True() : False()
@inline booltype(x::BoolsAsTypes)::BoolsAsTypes = x
