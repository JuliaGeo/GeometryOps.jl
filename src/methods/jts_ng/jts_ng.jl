# # JTS NG staging area
#
# This directory is the initial internal home for a Julia port of the JTS
# OverlayNG and RelateNG engines.  Keep the shared substrate small, and keep the
# RelateNG node model separate from the OverlayNG graph/label model.

include("common/topology.jl")
include("common/intersection_matrix.jl")
include("common/extraction.jl")
include("common/segment_primitives.jl")
include("common/algorithms.jl")

include("relateng/relateng.jl")
include("overlayng/overlayng.jl")
