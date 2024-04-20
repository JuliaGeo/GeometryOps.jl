# Predicates

Exact vs fast predicates

## Orient

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


## Incircle

