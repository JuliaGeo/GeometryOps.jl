module GeometryOpsFlexiJoinsExt

using GeometryOps
using FlexiJoins

import GeometryOps as GO, GeoInterface as GI
using SortTileRecursiveTree, Tables


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

FlexiJoins.prepare_for_join(::FlexiJoins.Mode.Tree, X, cond::FlexiJoins.ByPred{<: GO_DE9IM_FUNCS}) = (X, SortTileRecursiveTree.STRtree(cond.Rf(X)))
function FlexiJoins.findmatchix(::FlexiJoins.Mode.Tree, cond::FlexiJoins.ByPred{F}, ix_a, a, (B, tree)::Tuple, multi::typeof(identity)) where F <: GO_DE9IM_FUNCS
    idxs = SortTileRecursiveTree.query(tree, a)
    intersecting_idxs = filter(idxs) do idx
        cond.pred(a, cond.Rf(B)[idx])
    end
    return intersecting_idxs
end

# Finally, for completeness, we define the `swap_sides` function for those predicates which are defined as inversions.

FlexiJoins.swap_sides(::typeof(GO.contains)) = GO.within
FlexiJoins.swap_sides(::typeof(GO.within)) = GO.contains
FlexiJoins.swap_sides(::typeof(GO.coveredby)) = GO.covers
FlexiJoins.swap_sides(::typeof(GO.covers)) = GO.coveredby

# That's a wrap, folks!

end

