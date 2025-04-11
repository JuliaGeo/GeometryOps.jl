#=
# `BoolsAsTypes`

In `apply` and `applyreduce`, we pass `threading` and `calc_extent` as types, not simple boolean values.  

This is to help compilation - with a type to hold on to, it's easier for 
the compiler to separate threaded and non-threaded code paths.

Note that if we didn't include the parent abstract type, this would have been really 
type unstable, since the compiler couldn't tell what would be returned!

We had to add the type annotation on the `booltype(::Bool)` method for this reason as well.


!!! note Static.jl

    Static.jl is a package that provides a way to store and manipulate static values.
    But it creates a lot of invalidations since it breaks the assumption that operations
    like `<`, `>` and `==` can only return booleans.  So we don't use it here.

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

"""
    booltype(x)

Returns a [`BoolsAsTypes`](@ref) from `x`, whether it's a boolean or a BoolsAsTypes.
"""
function booltype end

@inline booltype(x::Bool)::BoolsAsTypes = x ? True() : False()
@inline booltype(x::BoolsAsTypes)::BoolsAsTypes = x

@inline istrue(x::True) = true
@inline istrue(x::False) = false
@inline istrue(x::Bool) = x
