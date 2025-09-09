# # TGGeometry extension
#=
[TGGeometry.jl](https://github.com/JuliaGeometry/TGGeometry.jl) is a Julia wrapper around the C library [`libtg`](https://github.com/tidwall/tg).
which has an innovative acceleration scheme for 2-D spatial predicates (intersects, contains, etc.)

This extension provides a GeometryOps interface to that using GeometryOps' [`TG`](@ref GO.TG) algorithm.

You can use any predicate like so:
```julia
import GeometryOps as GO
GO.intersects(GO.TG(), geom1, geom2)
```

or any other predicate in the list:

```@eval
using Markdown
import TGGeometry as TGG

Markdown.parse(
    join(["`$(f)`" for f in TGG.TG_PREDICATES], ", ")
)
```
=#
module GeometryOpsTGGeometryExt

using GeometryOps: TG
import GeometryOps as GO

using TGGeometry

# Literally loop over every name in TG_PREDICATES and eval in a function implementation.
# This is short and sweet, and completely static, so it doesn't run at runtime but rather
# at compile time.

# TODO: this could use some precompile statements, maybe.

for jl_fname in TGGeometry.TG_PREDICATES
    @eval GO.$jl_fname(::TG, geom1, geom2) = TGGeometry.$jl_fname(geom1, geom2)
end

end