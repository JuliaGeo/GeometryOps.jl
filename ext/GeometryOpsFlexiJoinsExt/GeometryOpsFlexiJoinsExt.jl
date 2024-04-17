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

FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.contains)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.within)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.intersects)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.disjoint)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.touches)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.crosses)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.overlaps)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.covers)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.coveredby)}, datas) = true
FlexiJoins.supports_mode(::FlexiJoins.Mode.NestedLoopFast, ::FlexiJoins.ByPred{typeof(GO.equals)}, datas) = true

# Next, just in case, we define the `swap_sides` function for those predicates which are defined as inversions.

FlexiJoins.swap_sides(::typeof(GO.contains)) = GO.within
FlexiJoins.swap_sides(::typeof(GO.within)) = GO.contains
FlexiJoins.swap_sides(::typeof(GO.coveredby)) = GO.covers
FlexiJoins.swap_sides(::typeof(GO.covers)) = GO.coveredby

# That's a wrap, folks!

end

