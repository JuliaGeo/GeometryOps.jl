
"""
    WrongManifoldError{InputManifold, DesiredManifold, Algorithm} <: Error

This error is thrown when an `Algorithm` is called with a manifold that it was not designed for.

It's mainly called for [`SingleManifoldAlgorithm`](@ref) types.
"""
struct WrongManifoldException{InputManifold, DesiredManifold, Algorithm} <: Base.Exception
    extra_text::String
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
    if !isempty(e.extra_text)
        print(io, "\n\n")
        print(io, e.extra_text)
    end
end