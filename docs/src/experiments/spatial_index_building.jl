using CairoMakie
import GeoInterface as GI, GeometryOps as GO
using NaturalEarth


using GeometryOps.SpatialTreeInterface
import GeometryOps.SpatialTreeInterface as STI
using SortTileRecursiveTree

function build_spatial_index_gif(geom, index_constructor, filename; plot_leaves = true, title = splitext(filename)[1], axis = (;), figure = (;), record = (;))
    fig = Figure(; figure...)
    ax = Axis(fig[1, 1]; title = title, axis...)

    # Create a spatial index
    index = index_constructor(geom)
    ext = STI.node_extent(index)
    limits!(ax, ext.X[1], ext.X[2], ext.Y[1], ext.Y[2])

    rects = Rect2f[Rect2f((NaN, NaN), (NaN, NaN))]
    colors = RGBAf[to_color(:transparent)]
    palette = Makie.wong_colors(0.7)

    plt = poly!(ax, rects; color = colors)

    to_rect2(extent) = Rect2f((extent.X[1], extent.Y[1]), (extent.X[2] - extent.X[1], extent.Y[2] - extent.Y[1]))

    function dive_in(io, plt, node, level)
        if STI.isleaf(node) && plot_leaves
            push!(rects, to_rect2(STI.node_extent(node)))
            push!(colors, palette[level])
        else
            for child in STI.getchild(node)
                dive_in(io, plt, child, level + 1)
            end
            push!(rects, to_rect2(STI.node_extent(node)))
            push!(colors, palette[level])
        end
        update!(plt, rects; color = colors)
        recordframe!(io)
        return
    end

    Makie.record(fig, filename; record...) do io
        empty!(rects)
        empty!(colors)
        dive_in(io, plt, index, 1)
    end

end 