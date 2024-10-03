#=
# Keyword docs

This file defines common keyword documentation, that can be spliced into docstrings.
=#

const THREADED_KEYWORD = "- `threaded`: `true` or `false`. Whether to use multithreading. Defaults to `false`."
const CRS_KEYWORD = "- `crs`: The CRS to attach to geometries. Defaults to `nothing`."
const CALC_EXTENT_KEYWORD = "- `calc_extent`: `true` or `false`. Whether to calculate the extent. Defaults to `false`."

const APPLY_KEYWORDS = """
$THREADED_KEYWORD
$CRS_KEYWORD
$CALC_EXTENT_KEYWORD
"""