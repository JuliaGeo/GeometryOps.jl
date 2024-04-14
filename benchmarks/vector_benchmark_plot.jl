using CairoMakie, MakieTeX

using CSV, DataFrames
using DataToolkit

path_to_makietex_datatoml = joinpath(dirname(dirname(@__DIR__)), "MakieTeX", "docs", "Data.toml")
data = DataToolkit.load(path_to_makietex_datatoml)


using DataToolkit, DataFrames, StatsBase
using CairoMakie, SwarmMakie #=beeswarm plots=#, Colors
using MakieTeX # for SVG icons

function svg_icon(name::String)
    if name == "go"
        icon = d"go-logo-solid::IO"
    else
        path = "svg/$name.svg"
        icon = get(d"file-icons::Dict{String,IO}", path, nothing)
    end
    if isnothing(icon)
        icon = get(d"file-icons-mfixx::Dict{String,IO}", path, nothing)
    end
    if isnothing(icon)
        icon = get(d"file-icons-devopicons::Dict{String,IO}", path, nothing)
    end
    isnothing(icon) && return missing
    return CachedSVG(read(seekstart(icon), String))
end

const colours_vibrant = range(LCHab(60,70,0), stop=LCHab(60,70,360), length=36)
const colours_dim     = range(LCHab(25,50,0), stop=LCHab(25,50,360), length=36)

const julia_logo = svg_icon("Julia")
const r_logo = svg_icon("R")
const python_logo = svg_icon("python")

marker_map = Dict(
    "geometryops" => julia_logo,
    # "gdal-jl" => julia_logo,
    "sf" => r_logo, 
    "terra" => r_logo, 
    "geos" => r_logo, 
    "s2" => r_logo,
    "geopandas" => python_logo,
)


color_map = Dict(
    # R packages
    "sf" => Makie.wong_colors()[1],
    "s2" => Makie.wong_colors()[5],
    "terra" => Makie.wong_colors()[6],
    "geos" => Makie.wong_colors()[4],
    # Python package
    "geopandas" => Makie.wong_colors()[2],
    # Julia package
    "geometryops" => Makie.wong_colors()[3],
)

path_to_vector_benchmark = "/Users/anshul/git/vector-benchmark"
timings_df = CSV.read(joinpath(path_to_vector_benchmark, "timings.csv"), DataFrame)
replace!(timings_df.package, "sf-project" => "sf", "sf-transform" => "sf")

# now plot

using SwarmMakie

using CategoricalArrays

task_ca = CategoricalArray(timings_df.task)

group_marker = [MarkerElement(; color = color_map[package], marker = marker_map[package], markersize = 12) for package in keys(marker_map)]
names_marker = collect(keys(marker_map))
lang_markers = ["R" => r_logo, "Python" => python_logo, "Julia" => julia_logo]
group_package = [MarkerElement(; marker, markersize = 12) for (lang, marker) in lang_markers]
names_package = first.(lang_markers)


f, a, p = beeswarm(
    task_ca.refs, timings_df.median;
    marker = getindex.((marker_map,), timings_df.package), 
    color = getindex.((color_map,), timings_df.package),
    markersize = 10,
    axis = (;
        xticks = (1:length(task_ca.pool.levels), task_ca.pool.levels),
        xlabel = "Task",
        ylabel = "Median time (s)",
        yscale = log10,
        title = "Benchmark vector operations",
        xgridvisible = false,
        xminorgridvisible = true,
        yminorgridvisible = true,
        yminorticks = IntervalsBetween(5),
        ygridcolor = RGBA{Float32}(0.0f0,0.0f0,0.0f0,0.05f0),
    )
)
leg = Legend(
    f[1, 2],
    [group_marker, group_package],
    [names_marker, names_package],
    ["Package", "Language"],
    tellheight = false,
    tellwidth = true,
    gridshalign = :left,
)
resize!(f, 650, 450)
a.spinewidth[] = 0.5
f