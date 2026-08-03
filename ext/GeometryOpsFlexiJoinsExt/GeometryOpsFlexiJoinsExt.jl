# # FlexiJoins extension
#=
This is an extension on FlexiJoins.jl that provides spatial tree acceleration to join 
arbitrary tables, using GeometryOps predicates.

The implementation is specialized on GO predicates.

See the [Spatial Joins tutorial](@ref tutorial-spatial-joins) for more information and for how to use this.
=#
module GeometryOpsFlexiJoinsExt

using GeometryOps
using FlexiJoins

import GeometryOps as GO, GeoInterface as GI
import GeometryOps.SpatialTreeInterface: spatialtree
using GeometryOps.SpatialTreeInterface: query, isspatialtree
using Tables


# This module defines the FlexiJoins APIs for GeometryOps' boolean comparison functions, taken from DE-9IM.

# First, we define the joining modes (Tree, NestedLoopFast) that the GO DE-9IM functions support.
const GO_DE9IM_FUNCS = Union{typeof(GO.contains), typeof(GO.within), typeof(GO.intersects), typeof(GO.disjoint), typeof(GO.touches), typeof(GO.crosses), typeof(GO.overlaps), typeof(GO.covers), typeof(GO.coveredby), typeof(GO.equals)}
# NestedLoopFast is the naive fallback method
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{F}, datas) where F <: GO_DE9IM_FUNCS = true
# This method allows you to cache a tree, which we do by using an STRtree.
# TODO: wrap GO predicate functions in a `TreeJoiner` struct or something, to indicate that we want to use trees,
# since they can be slower in some situations.
FlexiJoins.supports_mode(::FlexiJoins.Mode.Tree, ::FlexiJoins.ByPred{F}, datas) where F <: GO_DE9IM_FUNCS = true

# Nested loop support is simple, and needs no further support.  
# However, for trees, we need to define how the tree is prepared and how it is used.
# This is done by defining the `prepare_for_join` function to return an STRTree,
# and by defining the `findmatchix` function as querying that tree before checking
# intersections.

# In theory, one could extract the tree from e.g a GeoPackage or some future GeoDataFrame.

function spatialtree(X, selector)
    tree_or_geometries = selector(X)
    tree_or_geometries === nothing && return nothing
    ismissing(tree_or_geometries) && return nothing
    isspatialtree(tree_or_geometries) && return tree_or_geometries
    return spatialtree(tree_or_geometries)
end

FlexiJoins.prepare_for_join(::FlexiJoins.Mode.Tree, X, cond::FlexiJoins.ByPred{<: GO_DE9IM_FUNCS}) = (X, spatialtree(X, cond.Rf))
function FlexiJoins.findmatchix(::FlexiJoins.Mode.Tree, cond::FlexiJoins.ByPred{F}, ix_a, a, (B, tree)::Tuple, multi::typeof(identity)) where F<:GO_DE9IM_FUNCS
    # Implementation note:
    # here, `a` is a row, and `b` is the full table.
    # We extract the relevant columns using cond.Lf and cond.Rf.
    tree === nothing && return Int[]
    left_geom = cond.Lf(a)
    (isnothing(left_geom) || ismissing(left_geom) || GI.isempty(left_geom)) && return Int[]
    idxs = query(tree, left_geom)
    intersecting_idxs = filter!(idxs) do idx
        cond.pred(cond.Lf(a), cond.Rf(B[idx]))
    end
    return intersecting_idxs
end

# Finally, for completeness, we define the `swap_sides` function for those predicates which are defined as inversions.

FlexiJoins.swap_sides(::typeof(GO.contains)) = GO.within
FlexiJoins.swap_sides(::typeof(GO.within)) = GO.contains
FlexiJoins.swap_sides(::typeof(GO.coveredby)) = GO.covers
FlexiJoins.swap_sides(::typeof(GO.covers)) = GO.coveredby

function _geometry_column(table, input_index)
    geometry_columns = GI.geometrycolumns(table)
    if isnothing(geometry_columns) || isempty(geometry_columns)
        throw(ArgumentError("$input_index input does not declare a geometry column"))
    end
    geometry_column = first(geometry_columns)
    length(geometry_columns) > 1 && @warn "$input_index input declares multiple geometry columns $(repr(geometry_columns)); using the first."
    return geometry_column
end

function FlexiJoins.flexijoin((a, b)::NTuple{2,Any}, predicate::GO_DE9IM_FUNCS; kwargs...)
    left_col = _geometry_column(a, "First")
    right_col = _geometry_column(b, "Second")
    return FlexiJoins.flexijoin(
        (a, b),
        FlexiJoins.by_pred(left_col, predicate, right_col);
        kwargs...,
    )
end

# That's a wrap, folks!

end
