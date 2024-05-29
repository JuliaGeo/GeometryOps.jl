using CairoMakie, Chairmarks, BenchmarkTools
import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG
include("benchmark_plots.jl")

p25 = GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 8.0), (0.0, 8.0), (0.0, 0.0)], [(4.0, 0.5), (4.5, 0.5), (4.5, 3.5), (4.0, 3.5), (4.0, 0.5)], [(2.0, 4.0), (4.0, 4.0), (4.0, 6.0), (2.0, 6.0), (2.0, 4.0)]])
p26 = GI.Polygon([[(3.0, 1.0), (8.0, 1.0), (8.0, 7.0), (3.0, 7.0), (3.0, 5.0), (6.0, 5.0), (6.0, 3.0), (3.0, 3.0), (3.0, 1.0)], [(3.5, 5.5), (6.0, 5.5), (6.0, 6.5), (3.5, 6.5), (3.5, 5.5)],  [(5.5, 1.5), (5.5, 2.5), (3.5, 2.5), (3.5, 1.5), (5.5, 1.5)]])

suite = BenchmarkGroup(["title:Polygon intersection timing","subtitle:Single polygon, densified"])

for max_distance in exp10.(LinRange(-1, 1.5, 10))
    p25s = GO.segmentize(p25; max_distance)
    p26s = GO.segmentize(p26; max_distance)
    n_verts = GI.npoint(p25s)
    suite["GeometryOps"][n_verts] = @be GO.intersection($p25s, $p26s; target = $(GI.PolygonTrait()), fix_multipoly = $nothing)
    suite["LibGEOS"][n_verts] = @be LG.intersection($(GI.convert(LG, p25s)), $(GI.convert(LG, p26s)))
end


plot_trials(suite)