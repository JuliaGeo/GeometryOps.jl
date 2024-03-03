# # Performance benchmarks

# We'll benchmark GeometryOps.jl against LibGEOS, which is what most common geometry operations packages (across languages) seem to depend on.

# First, we'll load the desired packages:
import GeoInterface as GI, 
    GeometryBasics as GB,
    GeometryOps as GO, 
    LibGEOS as LG
import GeoInterface, GeometryBasics, GeometryOps, LibGEOS
using BenchmarkTools
using GeoInterface, GeometryBasics, GeoJSON # to generate and manipulate geometry
using CairoMakie # to visualize and understand what exactly we're doing
using DataInterpolations # to upscale and downscale geometry

# We set up a benchmark suite in order to understand exactly what will happen:
suite = BenchmarkGroup()

# In order to make this fair, we will each package's native representation as input to their benchmarks.
lg_and_gb(geometry) = (GI.convert(LibGEOS, geometry), GI.convert(GeometryBasics, geometry))

# ## Polygon benchmarks


circle_area_suite = BenchmarkGroup()
# Let's look at the simple case of a circle.
points = Point2f.((cos(θ) for θ in LinRange(0, 2π, 10000)), (sin(θ) for θ in LinRange(0, 2π, 10000)))

# We'll use this circle as a polygon for our benchmarks.
circle = GB.Polygon(points, [GB.decompose(Point2f, GB.Circle(Point2f(0.25, 0.25), 0.5))])
Makie.poly(circle; axis = (; aspect = DataAspect()))

closed_circle = ClosedRing()(circle)

lg_circle, go_circle = lg_and_gb(closed_circle)

# ### Area

# Let's start with the area of the circle.

@benchmark GO.area($go_circle)
@benchmark LG.area($lg_circle)

