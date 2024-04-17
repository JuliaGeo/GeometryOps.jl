module GeometryOpsFlexiJoinsExt

using GeometryOps
import GeometryOps as GO
using FlexiJoins

# This module defines the FlexiJoins APIs for GeometryOps' boolean comparison functions, taken from DE-9IM.

# For now, we only allow regular n^2 loops for the predicates, primarily because I'm not entirely sure how to do the "fast" mode for these.
# TODO: fix that and allow faster joins.

# First, we define that FlexiJoins supports the "NestedLoopFast" mode for the predicates.

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

