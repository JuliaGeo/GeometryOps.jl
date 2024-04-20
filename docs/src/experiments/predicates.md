# Predicates

Exact vs fast predicates

## Orient

```@example orient
using CairoMakie
import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG
import ExactPredicates
using MultiFloats

function orient_f64(p, q, r)
    return sign((GI.x(p) - GI.x(r))*(GI.y(q) - GI.y(r)) - (GI.y(p) - GI.y(r))*(GI.x(q) - GI.x(r)))
end

function orient_adaptive(p, q, r)
    px, py = Float64x2(GI.x(p)), Float64x2(GI.y(p))
    qx, qy = Float64x2(GI.x(q)), Float64x2(GI.y(q))
    rx, ry = Float64x2(GI.x(r)), Float64x2(GI.y(r))
    return sign((px - rx)*(qy - ry) - (py - ry)*(qx - rx))
end
# Create an interactive Makie dashboard which can show what is done here
fig = Figure()
axs = [Axis(fig[1, i]; aspect = DataAspect(), title) for (i, title) in enumerate(["Float64", "Adaptive", "Exact"])]
# function generate_heatmap()
```


### Dashboard
```julia
using WGLMakie
import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG
import ExactPredicates
using MultiFloats

function orient_f64(p, q, r)
    return sign((GI.x(p) - GI.x(r))*(GI.y(q) - GI.y(r)) - (GI.y(p) - GI.y(r))*(GI.x(q) - GI.x(r)))
end

function orient_adaptive(p, q, r)
    px, py = Float64x2(GI.x(p)), Float64x2(GI.y(p))
    qx, qy = Float64x2(GI.x(q)), Float64x2(GI.y(q))
    rx, ry = Float64x2(GI.x(r)), Float64x2(GI.y(r))
    return sign((px - rx)*(qy - ry) - (py - ry)*(qx - rx))
end
# Create an interactive Makie dashboard which can show what is done here
fig = Figure()
ax = Axis(fig[1, 1]; aspect = DataAspect())
sliders = SliderGrid(fig[2, 1],
        (label = L"w = 2^{-v} (zoom)", range = LinRange(40, 44, 100), startvalue = 42),
        (label = L"r = (x, y),~ x, y ∈ v + [0..w)", range = 0:0.01:3, startvalue = 0.95),
        (label = L"q = (k, k),~ k = v", range = LinRange(0, 30, 100), startvalue = 18),
        (label = L"p = (k, k),~ k = v", range = LinRange(0, 30, 100), startvalue = 16.8),
)
orient_funcs = [orient_f64, orient_adaptive, ExactPredicates.orient]
menu = Menu(fig[3, 1], options = zip(string.(orient_funcs), orient_funcs))
w_obs, r_obs, q_obs, p_obs = getproperty.(sliders.sliders, :value)
orient_obs = menu.selection

heatmap_size = @lift maximum(widths($(ax.scene.viewport)))*4

matrix_observable = lift(orient_obs, w_obs, r_obs, q_obs, p_obs, heatmap_size) do orient, w, r, q, p, heatmap_size
    return [orient((p, p), (q, q), (r+x, r+y)) for x in LinRange(0, 0+2.0^(-w), heatmap_size), y in LinRange(0, 0+2.0^(-w), heatmap_size)]
end
heatmap!(ax, matrix_observable; colormap = [:red, :green, :blue])
resize!(fig, 500, 700)
fig
```

### Testing robust vs regular predicates

```julia

import GeoInterface as GI, GeometryOps as GO
using MultiFloats
c1 = [[-28083.868447876892, -58059.13401805979], [-9833.052704767595, -48001.726711609794], [-16111.439295815226, -2.856614689791036e-11], [-76085.95770326033, -2.856614689791036e-11], [-28083.868447876892, -58059.13401805979]]
c2 = [[-53333.333333333336, 0.0], [0.0, 0.0], [0.0, -80000.0], [-60000.0, -80000.0], [-53333.333333333336, 0.0]]

p1 = GI.Polygon([c1])
p2 = GI.Polygon([c2])
GO.intersection(p1, p2; target = GI.PolygonTrait(), fix_multipoly = nothing)

p1_m, p2_m = GO.transform(x -> (Float64x2.(x)), [p1, p2])
GO.intersection(p1_m, p2_m; target = GI.PolygonTrait(), fix_multipoly = nothing)

p1 = GI.Polygon([[[-57725.80869813739, -52709.704377648755], [-53333.333333333336, 0.0], [-41878.01362848005, 0.0], [-36022.23699059147, -43787.61366192682], [-48268.44121252392, -52521.18593721105], [-57725.80869813739, -52709.704377648755]]])
p2 = GI.Polygon([[[-60000.0, 80000.0], [0.0, 80000.0], [0.0, 0.0], [-53333.33333333333, 0.0], [-50000.0, 40000.0], [-60000.0, 80000.0]]])
p1_m, p2_m = GO.transform(x -> (Float64x2.(x)), [p1, p2])
f, a, p__1 = poly(p1; label = "p1")
p__2 = poly!(a, p2; label = "p2")

GO.intersection(p1_m, p2_m; target = GI.PolygonTrait(), fix_multipoly = nothing)


```

## Incircle

