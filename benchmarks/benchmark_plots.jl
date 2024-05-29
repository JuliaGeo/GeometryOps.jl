using Printf, Statistics
using Makie, MakieThemes
using BenchmarkTools, Chairmarks

function compute_gridsize(numplts::Int, nr::Int, nc::Int)
    # figure out how many rows/columns we need
    if nr < 1
        if nc < 1
            nr = round(Int, sqrt(numplts))
            nc = ceil(Int, numplts / nr)
        else
            nr = ceil(Int, numplts / nc)
        end
    else
        nc = ceil(Int, numplts / nr)
    end
    nr, nc
end

function _prettytime(t::Real)
    if t < 1e-6
        value, units = t * 1e9, "ns"
    elseif t < 1e-3
        value, units = t * 1e6, "μs"
    elseif t < 1
        value, units = t * 1e3, "ms"
    else
        value, units = t * 1,   "s"
    end
    return string(@sprintf("%.1f", value), " ", units)
end

_prettytime(ts::AbstractArray{<: Real}) = _prettytime.(ts)


function results_to_numbers(result::BenchmarkGroup, postprocess_times = Statistics.median, postprocess_numbers = identity)
    # First, we extract the keys from the result.  
    # It's assumed that there is only one key per result, and that it's a number.
    numbers = identity.(collect(keys(result)))
    @assert numbers isa AbstractVector{<: Number} """
        Extra keys involved in the benchmark group!  
        Provide a pure group with only numerical keys.
        Got key types:
        $(unique(typeof.(numbers)))
        """
    # We sort the numbers, and then we use them to get the results
    sort!(numbers)
    result_objects = getindex.((result,), numbers)
    # Now, we get the times, and return their medians.
    postprocessed_objects = postprocess_times.(result_objects)
    return postprocess_numbers.(numbers), getproperty.(postprocessed_objects, :time)
end

# function plot_trials(results; title = )



"""
    plot_trials(result, title_str; theme = MakieThemes.bbc())::Figure
This function takes `result::BenchmarkTools.BenchmarkGroup` and plots the trials.  It returns the figure.
"""
plot_trials(results; kwargs...) = begin
    fig = Figure()
    plot_trials(fig[1, 1], results; kwargs...)
    fig
end

function plot_trials(
        gp::Makie.GridPosition,
        results;
        theme = merge(deepcopy(Makie.CURRENT_DEFAULT_THEME), MakieThemes.bbc()),
        legend_position = Makie.automatic, #(1, 1, TopRight()),
        legend_orientation = :horizontal,
        legend_halign = 1.0,
        legend_valign = -0.25,
    )

    xs, ys, labels = [], [], []
    for label in keys(results)
        current_result = results[label]
        if isempty(current_result)
            @warn "ResultSet with key $label is empty, skipping."
            continue
        end
        x, y = results_to_numbers(current_result)
        push!(xs, x)
        push!(ys, y)
        push!(labels, label)
    end

    tag_attrs = capture_tag_attrs(results.tags)

    lp = if legend_position isa Makie.Automatic 
        gp.layout[gp.span.rows, gp.span.cols, TopRight()] 
    elseif legend_position isa Tuple 
        gp.layout[legend_position...] 
    elseif legend_position isa Union{Makie.GridPosition, Makie.GridSubposition}
        legend_position
    else
        error()
    end

    ax = Makie.with_theme(theme) do
        ax = Axis(
            gp;
            tag_attrs.Axis...,
            xlabel = "Number of points", ylabel = "Time to calculate",
            xscale = log10, yscale = log10, ytickformat = _prettytime,
            xticksvisible = true, xticklabelsvisible = true,
            yticks = Makie.LogTicks(Makie.WilkinsonTicks(7; k_min = 4)),
            ygridwidth = 0.75,
        ) 
        plots = [scatterlines!(ax, x, y; label = label) for (x, y, label) in zip(xs, ys, labels)]
        setproperty!.(getindex.(getproperty.(plots, :plots), 1), :alpha, 0.1)
        leg = Legend(
            lp, ax; 
            tellwidth = legend_position isa Union{Tuple, Makie.Automatic} && (legend_position isa Makie.Automatic || length(legend_position) != 3) && legend_orientation == :vertical, 
            tellheight = legend_position isa Union{Tuple, Makie.Automatic} && (legend_position isa Makie.Automatic || length(legend_position) != 3) && legend_orientation == :horizontal, 
            halign = legend_halign, 
            valign = legend_valign, 
            orientation = legend_orientation
        )
        ax.xtickcolor[] = ax.xgridcolor[]
        ax
    end

    return ax
end

const _tag_includelist = ["title", "subtitle"]

function capture_tag_attrs(tags)
    attr_dict = Attributes()
    axis = attr_dict.Axis = Attributes()
    for tag in tags
        for possibility in sort(_tag_includelist; by = length)
            if startswith(tag, possibility)
                axis[Symbol(possibility)] = tag[(length(possibility) + 2):end]
                break
            end
        end
    end
    return attr_dict
end

function decompose_benchmarksuite_to_2ndlevel(result)
    # here, `result` is a BenchmarkGroup.
    
end