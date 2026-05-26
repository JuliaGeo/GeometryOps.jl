# RelateNG-specific code goes here.
#
# Keep RelateNG node sections, node edges, topology predicates, and point
# locators out of `common/` unless their semantics are genuinely shared with
# OverlayNG.

include("point_locator.jl")
include("relate_geometry.jl")
include("topology_predicates.jl")
