using Printf, Statistics, MakieThemes

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


function results_to_numbers(result::BenchmarkTools.BenchmarkGroup, postprocess_times = Statistics.median, postprocess_numbers = identity)
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
function plot_trials(
        results, 
        title_str_part; 
        theme = merge(Makie.CURRENT_DEFAULT_THEME, MakieThemes.bbc()),
        legend_position = (1, 1, TopRight()),
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

    return Makie.with_theme(theme) do
        fig = Figure()
        ax = Axis(
            fig[1, 1];
            title = title_str_part, subtitle = "Tested on a regular circle",
            xlabel = "Number of points", ylabel = "Time to calculate",
            xscale = log10, yscale = log10, ytickformat = _prettytime,
        ) 
        plots = [scatterlines!(ax, x, y; label = label) for (x, y, label) in zip(xs, ys, labels)]
        setproperty!.(getindex.(getproperty.(plots, :plots), 1), :alpha, 0.1)
        leg = Legend(
            fig[legend_position...], ax; 
            tellwidth = length(legend_position) != 3 && legend_orientation == :vertical, 
            tellheight = length(legend_position) != 3 && legend_orientation == :horizontal, 
            halign = legend_halign, 
            valign = legend_valign, 
            orientation = legend_orientation
        )
        ax.xticksvisible[] = true
        ax.xtickcolor[] = ax.xgridcolor[]
        ax.xticklabelsvisible[] = true
        ax.yticks[] = Makie.LogTicks(Makie.WilkinsonTicks(7; k_min = 4))
        ax.ygridwidth[] = 0.75
        return fig
    end
end