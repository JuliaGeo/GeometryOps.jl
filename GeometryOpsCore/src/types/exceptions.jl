#=
# Errors and exceptions

We create a few custom exception types in this file,
that have nice show methods that we can use for certain errors.

This makes it substantially easier to catch specific kinds of errors and show them.
For example, we can catch `WrongManifoldException` and show a nice error message,
and error hinters can be specialized to that as well.

```@docs; canonical=false
GeometryOpsCore.WrongManifoldException
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