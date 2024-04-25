#=
# Geometry providers

This file benchmarks GeometryOps methods on every GeoInterface.jl implementation we can find, in order to test:
a. genericness, i.e., does GeometryOps work correctly with all GeoInterface.jl implementations?
b. performance, i.e., how does GeometryOps compare to the native implementation?
c. performance issues in the packages' implementations of GeoInterface
=#

# First, we import the providers:
using ArchGDAL, LibGEOS, Shapefile, GeoJSON, WellKnownGeometry, GeometryBasics, GeoInterface, GeoFormatTypes
PROVIDERS = (ArchGDAL, LibGEOS, GeometryBasics, GI.Wrappers)
# Now, we import GeoInterface and GeometryOps,
import GeometryOps as GO, GeoInterface as GI
# Finally, we import some utility benchmarking, plotting and data munging packages!
using BenchmarkTools, Chairmarks, CairoMakie, MakieThemes, DataFrames, Proj
using CoordinateTransformations, Rotations


# Polylabel.jl is a package that finds the "pole of inaccessibility" of a polygon,
# i.e., the point within it that is furthest away from its boundaries.  

# It depends on GeometryOps, but in this instance, we'll grab some of its test geometries
# to use.
import Polylabel

# TODO: the reason we change to LibGEOS intermediately here is so that the 
# linear rings of the WKG polygons are interpreted correctly.  Unfortunately 
# that doesn't work when read, which there's an issue up for.
water1 = GeoFormatTypes.WellKnownText(GeoFormatTypes.Geom(), readchomp(joinpath(dirname(dirname(pathof(Polylabel))), "test", "data", "water1.wkt")) |> String) |> x -> GI.convert(LibGEOS, x) |> GO.tuples
water2 = GeoFormatTypes.WellKnownText(GeoFormatTypes.Geom(), readchomp(joinpath(dirname(dirname(pathof(Polylabel))), "test", "data", "water2.wkt")) |> String) |> x -> GI.convert(LibGEOS, x) |> GO.tuples
# To fix these polygons is a complicated task, and even then LibGEOS gets it wrong:
# water1 |> x -> LibGEOS.makeValid(GI.convert(LibGEOS, x)) |> GI.getgeom |> collect |> x -> filter(y -> GI.trait(y) isa Union{GI.PolygonTrait, GI.MultiPolygonTrait}, x) |> first |> GO.tuples # hide

f, a, p = poly(water1; axis = (; title = "water1")); poly(f[1, 2], water2; axis = (; title = "water2")); f
# Now, we rotate the `water1` polygon about its centroid, so we can use it to 
# test the time it takes to intersect complex polygons:
water1r = GO.transform(
    Translation(GO.centroid(water1)) ∘ LinearMap(Makie.rotmatrix2d(π/2)) ∘ Translation((-).(GO.centroid(water1))), 
    water1
)
f, a, p = poly(water1; label = "Original")
poly!(water1r; label = "Rotated")
axislegend(a)
f
# WARNING: does not work
@b GO.union($(water1), $(water1r); target = GI.PolygonTrait()) seconds=3
@b LibGEOS.union($(GI.convert(LibGEOS, water1)), $(GI.convert(LibGEOS, water1r))) seconds=3
@b ArchGDAL.union($(GI.convert(ArchGDAL, water1)), $(GI.convert(ArchGDAL, water1r))) seconds=3

poly(GO.union(w1g, w1rg; target = GI.PolygonTrait()))

GI.getgeom(water1, 3) |> GI.trait

# We can benchmark each provider and see if any of them have glaring issues.

water1_centroid_suite = BenchmarkGroup()

for provider in PROVIDERS
    @info "Benchmarking $provider"
    geom = GI.convert(provider, water1)
    water1_centroid_suite[string(provider)] = @be GO.centroid($geom) seconds=3
end


# ## Tables.jl performance in `apply`
#=
This code checks how Tables.jl performs when using `apply`.
We use two sources for this: `Shapefile.jl` and `DataFrames.jl`.
More will be coming in the future!
=#
shp_file = "/Users/anshul/Downloads/ne_10m_admin_0_countries (1)/ne_10m_admin_0_countries.shp"
table = Shapefile.Table(shp_file)
go_df = DataFrame(table)
go_df.geometry = GO.tuples(go_df.geometry);

table_suite = BenchmarkGroup()


ll2moll = Proj.Transformation("+proj=longlat +datum=WGS84", "+proj=moll")

# First, we try reprojecting the geometries using Proj,
reproject_suite = table_suite["reproject"] = BenchmarkGroup(["title:Reproject", "subtitle:All country borders from Natural Earth, 1:10m res."])

reproject_suite["Shapefile.Table"] = @be GO.reproject($table, $ll2moll) seconds=3
reproject_suite["DataFrame (Shapefile)"] = @be GO.reproject($(DataFrame(table)), $ll2moll) seconds=3
reproject_suite["DataFrame (GO)"] = @be GO.reproject($(go_df), $ll2moll) seconds=3
reproject_suite["Shapefile geoms"] = @be GO.reproject($(table.geometry), $ll2moll) seconds=3
reproject_suite["GeometryOps geoms"] = @be GO.reproject($(GO.tuples(table.geometry)), $ll2moll) seconds=3

# then transforming, just to see the difference in runtime
# between calling out to C vs pure Julia,
function _scaleby5(x)
    return x .* 5
end

transform_suite = table_suite["transform"] = BenchmarkGroup(["title:Transform", "subtitle:All country borders from Natural Earth, 1:10m res."])
transform_suite["Shapefile.Table"] = @be GO.transform($_scaleby5, $table) seconds=3
transform_suite["DataFrame (Shapefile)"] = @be GO.transform($_scaleby5, $(DataFrame(table))) seconds=3
transform_suite["DataFrame (GO)"] = @be GO.transform($_scaleby5, $(go_df)) seconds=3
transform_suite["Shapefile geoms"] = @be GO.transform($_scaleby5, $(table.geometry)) seconds=3
transform_suite["GeometryOps geoms"] = @be GO.transform($_scaleby5, $(GO.tuples(table.geometry))) seconds=3

# and finally, calling `applyreduce` to find the area of each
# polygon.
area_suite = table_suite["area"] = BenchmarkGroup(["title:Area", "subtitle:All country borders from Natural Earth, 1:10m res."])

area_suite["Shapefile.Table"] = @be GO.area($(table)) seconds=3
area_suite["DataFrame (Shapefile)"] = @be GO.area($(DataFrame(table))) seconds=3
area_suite["DataFrame (GO)"] = @be GO.area($(go_df)) seconds=3
area_suite["Shapefile geoms"] = @be GO.area($(table.geometry)) seconds=3
area_suite["GeometryOps geoms"] = @be GO.area($(GO.tuples(table.geometry))) seconds=3

ts = getproperty.(area_suite["Shapefile.Table"].samples, :time)
boxplot(ones(length(ts)), ts)
violin(ones(length(ts)), ts; npoints = 3500, axis = (; yscale = log10,))


# ## Plotting
function Makie.convert_arguments(::Makie.PointBased, xs, bs::AbstractVector{<: Chairmarks.Benchmark})
    ts = getproperty.(Statistics.mean.(bs), :time)
    return (xs, ts)
end

function Makie.convert_arguments(::Makie.PointBased, bs::AbstractVector{<: Chairmarks.Benchmark})
    ts = getproperty.(Statistics.mean.(bs), :time)
    return (1:length(bs), ts)
end

function Makie.convert_arguments(::Makie.SampleBased, b::Chairmarks.Benchmark)
    ts = getproperty.(b.samples, :time)
    return (ones(length(ts)), ts)
end

function Makie.convert_arguments(::Makie.SampleBased, n::Number, b::Chairmarks.Benchmark)
    ts = getproperty.(b.samples, :time)
    return (fill(n, length(ts)), ts)
end

function Makie.convert_arguments(::Makie.SampleBased, labels::AbstractVector{<: AbstractString}, bs::AbstractVector{<: Chairmarks.Benchmark})
    ts = map(b -> getproperty.(b.samples, :time), bs)
    labels = 
    return flatten
end

function Makie.convert_arguments(::Type{Makie.Errorbars}, xs, bs::AbstractVector{<: Chairmarks.Benchmark})
    ts = map(b -> getproperty.(b.samples, :time), bs)
    means = map(Statistics.mean, ts)
    stds = map(Statistics.std, ts)
    return (xs, ts)
end

ks = keys(area_suite) |> collect .|> identity

bs = getindex.((area_suite,), ks)
b_lengths = length.(getproperty.(bs, :samples))
b_timing_flattened = collect(Iterators.flatten(Iterators.map(b -> getproperty.(b.samples, :time), bs)))
k_strings = Iterators.flatten((fill(k, bl) for (k, bl) in zip(ks, b_lengths))) |> collect

f = Figure()
ax = Axis(f[1, 1];
    convert_dim_1=Makie.CategoricalConversion(; sortby=nothing),
)
violin!(ax, k_strings, b_timing_flattened .|> log10)
f
ax.yscale = log10
ax.xticklabelrotation = π/12
f


bs = values(area_suite) |> collect .|> identity
labels = ["ST", "DS", "DG", "SG", "GG"]


using AlgebraOfGraphics

boxplot(b1)
boxplot!.(1:5, values(area_suite) |> collect .|> identity)
Makie.current_figure()
Makie.current_axis().yscale = log10

data((; x = labels, y = bs)) * mapping(:y => verbatim, :x, :y) * visual(BoxPlot) |> draw
