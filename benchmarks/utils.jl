using Printf, Statistics

function _prettytime(t::Real)
    if t < 1e3
        value, units = t, "ns"
    elseif t < 1e6
        value, units = t / 1e3, "μs"
    elseif t < 1e9
        value, units = t / 1e6, "ms"
    else
        value, units = t / 1e9, "s"
    end
    return string(@sprintf("%.1f", value), " ", units)
end

_prettytime(ts::AbstractArray{<: Real}) = _prettytime.(ts)

"""
    results_to_numbers(result::BenchmarkTools.BenchmarkGroup, postprocess = Statistics.median)

This function takes a `BenchmarkTools.BenchmarkGroup` and returns the numbers and the results of the group.  
    
The `postprocess` function is applied to the times vector, and by default it's `Statistics.median`.

"""
function results_to_numbers(result::BenchmarkTools.BenchmarkGroup, postprocess = Statistics.median)
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
    time_vectors = getproperty.(result_objects, :times)
    return numbers, postprocess.(time_vectors)
end

"""
    plot_trials(result, title_str; theme = MakieThemes.bbc())::Figure

This function takes `result::BenchmarkTools.BenchmarkGroup` and plots the trials.  It returns the figure.
"""
plot_trials(result::BenchmarkTools.BenchmarkGroup, title_str; theme = MakieThemes.bbc()) = plot_trials(result["LibGEOS"], result["GeometryBasics"], title_str)

function plot_trials(libgeos_results, geometryops_results, title_str_part; theme = MakieThemes.bbc())

    x_libgeos, y_libgeos = results_to_numbers(libgeos_results)
    x_geometryops, y_geometryops = results_to_numbers(geometryops_results)

    return Makie.with_theme(theme) do
        f, a, p1 = scatterlines(x_libgeos, y_libgeos; label = "LibGEOS", axis = (; xscale = log10, yscale = log10, ytickformat = _prettytime))
        p2 = scatterlines!(a, x_geometryops, y_geometryops; label = "GeometryOps")
        leg = Legend(f[1, 1, TopRight()], a; tellwidth = false, tellheight = false, halign = 1.0, valign = -0.25, orientation = :horizontal)
        for plot in (p1, p2)
            plot.plots[1].alpha[] = 0.1
        end 
        a.title = "$title_str_part"
        a.xlabel = "Number of points"
        a.ylabel = "Time to calculate"
        a.xticksvisible[] = true
        a.xtickcolor[] = a.xgridcolor[]
        a.xticklabelsvisible[] = true
        a.yticks[] = Makie.LogTicks(Makie.WilkinsonTicks(7; k_min = 4))
        a.ygridwidth[] = 0.75
        a.subtitle = "Tested on a regular circle"
        return f
    end
end