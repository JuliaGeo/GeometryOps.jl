#=
# A random walk through [GeometryOps.jl](https://juliageo.org/GeometryOps.jl)

In this tutorial, we'll take a random walk through GeometryOps and its capabilities, just to show off what it can do.


=#
import GeometryOps as GO
## ecosystem packages we'll need
import GeoInterface as GI
using CairoMakie # for plotting
#

#=
## The `apply` interface

- my_coord_op
- my_linestring_op
=# 


#=
## The `applyreduce` interface
This one is arguably more useful for daily tasks.

- my_centroid on a linestring/ring level
- 

=#

#=
## The `fix` interface
- Choose your fixes
- How to make a new fix (antimeridian cutting)
=#


#=
## LibGEOS extension
> If you can't do it yourself, then use something else.  
TODO: chatgpt this quote

=#
import LibGEOS # we will never actually call LibGEOS here

GO.buffer(poly, 1) |> Makie.poly
