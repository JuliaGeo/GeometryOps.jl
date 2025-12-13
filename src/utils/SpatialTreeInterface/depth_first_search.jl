#=
# Depth-first search

```@meta
CollapsedDocStrings = true
```

```@docs; canonical=false
depth_first_search
```

This file houses a function that performs a depth-first search over a spatial tree, filtering
by `predicate` on the extents of the nodes, and calling `f` on every leaf node that matches 
`predicate`.

Here, `predicate` is a one-argument function that returns a boolean, and `f` is a function 
that takes a single argument, the index of a geometry as stored in the leaf node.

The `depth_first_search` function is generic to anything that implements the SpatialTreeInterface.

## Example

Here's an animated example of what `depth_first_search` is doing:


```julia
using Makie
using Extents, GeometryOps.SpatialTreeInterface, SortTileRecursiveTree
import GeoInterface as GI

# Construct a tree of extents (in this case, just a grid)
xs, ys = 1.0:0.1:10.0, 1.0:0.1:10.0
rects_vec = vec([GI.Polygon([[(x, y), (x+0.1, y), (x+0.1, y+0.1), (x, y+0.1), (x, y)]]) for x in xs, y in ys])
ext2rect(ext) = Rect2f((ext.X[1], ext.Y[1]), (ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1]))

# First, create a figure that holds some plots:
f, a, all_rects_plot = lines(rects_vec; linewidth = 0.5, axis = (; aspect = DataAspect()));
Legend(f[1, 2], [Makie.PolyElement(; color = :transparent, strokecolor = :blue, strokewidth = 0.5), Makie.PolyElement(; color = :green), Makie.PolyElement(; color = :red)], ["Grid", "Branch node", "Leaf node"])

target_extent = Extents.Extent(X = (5, 7), Y = (5, 7))
target_rect_plot = poly!(a, ext2rect(target_extent); color = (:blue, 0.5))

current_rect = Observable(ext2rect(GI.extent(rects_vec[1])))
current_color = Observable(:red)

current_rect_plot = poly!(a, current_rect; color = current_color, alpha = 0.3)


record(f, "./depth_first_search.mp4"; framerate = 1) do io

    function leaf_f(i)
        current_rect[] = ext2rect(GI.extent(rects_vec[i]))
        current_color[] = :red
        recordframe!(io)
        return nothing
    end

    function pred_f(ext, target)
        current_rect[] = ext2rect(ext)
        current_color[] = :green
        recordframe!(io)
        result = Extents.intersects(ext, target)
        if result
            current_color[] = :forestgreen
            recordframe!(io)
        end
        return result
    end

    SpatialTreeInterface.depth_first_search(leaf_f, Base.Fix2(pred_f, target_extent), extents_tree)
end

````

=#


"""
    depth_first_search(f, predicate, tree)

Call `f(i)` for each index `i` in the tree that satisfies `predicate(extent(i))`.

This is generic to anything that implements the SpatialTreeInterface, particularly the methods
[`isleaf`](@ref), [`getchild`](@ref), and [`child_extents`](@ref).

## Example

```jldoctest
using Extents, GeometryOps.SpatialTreeInterface, SortTileRecursiveTree
# Construct a tree of extents (in this case, just a grid)
xs, ys = 1:10, 1:10
extents_vec = vec([Extents.Extent(X = (x, x+1), Y = (y, y+1)) for x in xs, y in ys])

# Construct a tree of extents - this is an STRtree,
# but works on any SpatialTreeInterface-compatible tree.
extents_tree = STRtree(extents_vec)

# Count the number of extents that intersect the extent (1,1) x (2,2)
target_extent = Extents.Extent(X = (1,2), Y = (1,2))
count = 0
SpatialTreeInterface.depth_first_search(
    (i::Int -> global count += 1), # `f`
    Base.Fix1(Extents.intersects, target_extent), 
    extents_tree
)
count
# output
4
```
"""
function depth_first_search(f::F, predicate::P, node::N) where {F, P, N}
    if isleaf(node)
        for (i, leaf_geometry_extent) in child_indices_extents(node)
            if predicate(leaf_geometry_extent)
                @controlflow f(i)
            end
        end
    else
        for child in getchild(node)
            if predicate(node_extent(child))
                @controlflow depth_first_search(f, predicate, child)
            end
        end
    end
end
function depth_first_search(predicate, node)
    a = Int[]
    depth_first_search(Base.Fix1(push!, a), predicate, node)
    return a
end
