#=
# Errors and exceptions

We create a few custom exception types in this file,
that have nice show methods that we can use for certain errors.

This makes it substantially easier to catch specific kinds of errors and show them.
For example, we can catch `WrongManifoldException` and show a nice error message,
and error hinters can be specialized to that as well.

We also have a custom error type for missing keywords in an algorithm,
which could eventually be extended to have typo detection, et cetera.

```@docs; canonical=false
WrongManifoldException
MissingKeywordInAlgorithmException
```

```@meta
CollapsedDocStrings = true
```
=#

export WrongManifoldException

"""
    WrongManifoldException{InputManifold, DesiredManifold, Algorithm} <: Exception

This error is thrown when an `Algorithm` is called with a manifold that it was not designed for.

It's mainly thrown when constructing `SingleManifoldAlgorithm` types.
"""
struct WrongManifoldException{InputManifold, DesiredManifold, Algorithm} <: Base.Exception
    description::String
end

WrongManifoldException{I, D, A}() where {I, D, A} = WrongManifoldException{I, D, A}("")

function Base.showerror(io::IO, e::WrongManifoldException{I,D,A}) where {I,D,A}
    print(io, "Algorithm ")
    printstyled(io, A; bold = true, color = :green) 
    print(io, " is only compatible with manifold ")
    printstyled(io, D; bold = true, color = :blue)
    print(io, ",\n but it was called with manifold ")
    printstyled(io, I; bold = true, color = :red)
    print(io, ".")

    println(io, """
    \n
    To fix this issue, you can specify the manifold explicitly, 
    e.g. `$A($D(); kwargs...)`, when constructing the algorithm.
    """)
    if !isempty(e.description)
        print(io, "\n\n")
        print(io, e.description)
    end
end


"""
    MissingKeywordInAlgorithmException{Alg, F} <: Exception

An error type which is thrown when a keyword argument is missing from an algorithm.

The `alg` argument is the algorithm struct, and the `keyword` argument is the keyword
that was missing.

This error message is used in the [`enforce`](@ref) method.

## Usage

This is of course not how you would actually use this error type, but it is
how you construct and throw it.

```julia
throw(MissingKeywordInAlgorithmException(GEOS(; tokl = 1.0), my_function, :tol))
```

Real world usage will often look like this:

```julia
function my_function(alg::CLibraryPlanarAlgorithm, args...)
    mykwarg = enforce(alg, :mykwarg, my_function) # this will throw an error if :mykwarg is not present in alg
end
```
"""
struct MissingKeywordInAlgorithmException{Alg, F} <: Base.Exception
    alg::Alg
    f::F
    keyword::Symbol
end

_name_of(x::Any) = nameof(typeof(x))
_name_of(x::Function) = nameof(x)
_name_of(x::Type) = nameof(x)
_name_of(s::String) = s

# This is just the generic dispatch, different algorithms can choose to dispatch
# on their own types to provide more specific or interesting error messages.
function Base.showerror(io::IO, e::MissingKeywordInAlgorithmException)
    algorithm_name = _name_of(typeof(e.alg))
    function_name = _name_of(e.f)
    print(io, "The ")
    printstyled(io, e.keyword; color = :red)
    print(io, " parameter is required for the ")
    printstyled(io, algorithm_name; bold = true)
    println(io, " algorithm in `$(function_name)`,")
    println(io, "but it was not provided.")
    println(io)
    println(io, "Provide it to the algorithm at construction time, like so:")
    println(io, "`$(algorithm_name)(; $(e.keyword) = ...)`")
    println(io, "and pass that as the algorithm to `$(function_name)`, usually the first argument.")
end