import GeometryOps as GO
# ecosystem packages we'll need
import GeoInterface as GI
using CairoMakie # for plotting


#=
## LibGEOS extension
> If you can't do it yourself, then use something else.  
TODO: chatgpt this quote

=#
import LibGEOS # we will never actually call LibGEOS here

GO.buffer(poly, 1) |> Makie.poly
