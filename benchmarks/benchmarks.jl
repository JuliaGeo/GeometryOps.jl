# # Performance benchmarks

# We'll benchmark GeometryOps.jl against LibGEOS, which is what most common geometry operations packages (across languages) seem to depend on.

# First, we'll load the desired packages:
import GeoInterface as GI, 
    GeometryBasics as GB,
    GeometryOps as GO, 
    LibGEOS as LG
import GeoInterface, GeometryBasics, GeometryOps, LibGEOS
using BenchmarkTools, Statistics
using GeoJSON # to generate and manipulate geometry
using CairoMakie, MakieThemes, GeoInterfaceMakie # to visualize and understand what exactly we're doing
using DataInterpolations # to upscale and downscale geometry

GeoInterfaceMakie.@enable GeoJSON.AbstractGeometry
GeoInterfaceMakie.@enable LibGEOS.AbstractGeometry
GeoInterfaceMakie.@enable GeoInterface.Wrappers.WrapperGeometry


# We include some basic plotting utilities here!
include(joinpath(@__DIR__, "utils.jl"))

# We set up a benchmark suite in order to understand exactly what will happen:
suite = BenchmarkGroup()

# In order to make this fair, we will each package's native representation as input to their benchmarks.
lg_and_go(geometry) = (GI.convert(LibGEOS, geometry), GO.tuples(geometry))

# # Polygon benchmarks


# Let's look at the simple case of a circle.
points = Point2f.((cos(θ) for θ in LinRange(0, 2π, 10000)), (sin(θ) for θ in LinRange(0, 2π, 10000)))

# We'll use this circle as a polygon for our benchmarks.
circle = GI.Wrappers.Polygon([points, GB.decompose(Point2f, GB.Circle(Point2f(0.25, 0.25), 0.5))])
closed_circle = GO.ClosedRing()(GO.tuples(circle))
Makie.poly(circle; axis = (; aspect = DataAspect()))
# Now, we materialize our LibGEOS circles;
lg_circle, go_circle = lg_and_go(closed_circle)

# ## Area

# Let's start with the area of the circle.
circle_area_suite = BenchmarkGroup()

# We compute the area of the circle at different resolutions!
n_points_values = [10, 100, 1000, 10000, 100000]
for n_points in n_points_values
    circle = GI.Wrappers.Polygon([tuple.((cos(θ) for θ in LinRange(0, 2π, n_points)), (sin(θ) for θ in LinRange(0, 2π, n_points)))])
    closed_circle = GO.ClosedRing()(circle)
    lg_circle, go_circle = lg_and_go(closed_circle)
    circle_area_suite["GeometryOps"][n_points] = @benchmarkable GO.area($go_circle)
    circle_area_suite["LibGEOS"][n_points]     = @benchmarkable LG.area($lg_circle)
end

BenchmarkTools.tune!(circle_area_suite)
circle_area_result = BenchmarkTools.run(circle_area_suite)

# We now have the benchmark results, and we can visualize them.

plot_results(circle_area_result, "Area")

# ## Difference, intersection, union

circle_suite = BenchmarkGroup()
circle_difference_suite = circle_suite["difference"]
circle_intersection_suite = circle_suite["intersection"]
circle_union_suite = circle_suite["union"]

n_points_values = round.(Int, exp10.(LinRange(1, 4, 10)))
for n_points in n_points_values
    circle = GI.Wrappers.Polygon([tuple.((cos(θ) for θ in LinRange(0, 2π, n_points)), (sin(θ) for θ in LinRange(0, 2π, n_points)))])
    closed_circle = GO.ClosedRing()(circle)

    lg_circle_right, go_circle_right = lg_and_go(closed_circle)

    circle_left = GO.apply(GI.PointTrait, closed_circle) do point
        x, y = GI.x(point), GI.y(point)
        return (x+0.6, y)
    end
    lg_circle_left, go_circle_left = lg_and_go(circle_left)
    circle_difference_suite["GeometryOps"][n_points] = @benchmarkable GO.difference($go_circle_left, $go_circle_right; target = GI.PolygonTrait)
    circle_difference_suite["LibGEOS"][n_points]     = @benchmarkable LG.difference($lg_circle_left, $lg_circle_right)
    circle_intersection_suite["GeometryOps"][n_points] = @benchmarkable GO.intersection($go_circle_left, $go_circle_right; target = GI.PolygonTrait)
    circle_intersection_suite["LibGEOS"][n_points]     = @benchmarkable LG.intersection($lg_circle_left, $lg_circle_right)
    circle_union_suite["GeometryOps"][n_points] = @benchmarkable GO.union($go_circle_left, $go_circle_right; target = GI.PolygonTrait)
    circle_union_suite["LibGEOS"][n_points]     = @benchmarkable LG.union($lg_circle_left, $lg_circle_right)
end

BenchmarkTools.tune!(circle_suite)
@time circle_result = BenchmarkTools.run(circle_suite; seconds = 3)

# Now, we plot!

# ### Difference
plot_trials(circle_result["difference"], "Difference")

# ### Intersection
plot_trials(circle_result["intersection"], "Intersection")

# ### Union
plot_trials(circle_result["union"], "Union")


# ## Good old USA
fc = GeoJSON.read(read(download("https://rawcdn.githack.com/nvkelso/natural-earth-vector/ca96624a56bd078437bca8184e78163e5039ad19/geojson/ne_10m_admin_0_countries.geojson")))
usa_multipoly = fc.geometry[findfirst(==("United States of America"), fc.NAME)]
areas = [GO.area(p) for p in GI.getgeom(usa_multipoly)]
usa_poly = GI.getgeom(usa_multipoly, findmax(areas)[2])
center_of_the_world = GO.centroid(usa_poly)
usa_poly_reflected = GO.apply(GI.PointTrait, usa_poly) do point
    x, y = GI.x(point), GI.y(point)
    return (-(x - GI.x(center_of_the_world)) + GI.x(center_of_the_world), y)
end

f, a, p = poly(usa_poly; color = Makie.wong_colors(0.5)[1], label = "Straight", axis = (; title = "Good old U.S.A.", aspect = DataAspect()))
poly!(usa_poly_reflected; color = Makie.wong_colors(0.5)[2], label = "Reversed")
Legend(f[2, 1], a; valign = 0, orientation = :horizontal)
f

usa_o_lg, usa_o_go = lg_and_go(usa_poly)
usa_r_lg, usa_r_go = lg_and_go(usa_poly_reflected)

# First, we'll test union:
printstyled("LibGEOS"; color = :red, bold = true)
println()
@benchmark LG.union($usa_o_lg, $usa_r_lg)
printstyled("GeometryOps"; color = :blue, bold = true)
println()
@benchmark GO.union($usa_o_go, $usa_r_go; target = GI.PolygonTrait)

# Next, intersection:
printstyled("LibGEOS"; color = :red, bold = true)
println()
@benchmark LG.intersection($usa_o_lg, $usa_r_lg)
printstyled("GeometryOps"; color = :blue, bold = true)
println()
@benchmark GO.intersection($usa_o_go, $usa_r_go; target = GI.PolygonTrait)

# and finally the difference:
printstyled("LibGEOS"; color = :red, bold = true)
println()
lg_diff = LG.difference(usa_o_lg, usa_r_lg)
printstyled("GeometryOps"; color = :blue, bold = true)
println()
go_diff = GO.difference(usa_o_go, usa_r_go; target = GI.PolygonTrait)

# You can see clearly that GeometryOps is currently losing out to LibGEOS.  Our algorithms aren't optimized for large polygons and we're paying the price for that.

# It's heartening that the polygon complexity isn't making too much of a difference; the difference in performance is mostly due to the number of vertices, as we can see from the circle benchmarks as well.