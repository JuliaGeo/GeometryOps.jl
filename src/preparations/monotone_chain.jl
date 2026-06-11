#=
# Monotone chain

A monotone chain is a continuous list of edges whose slopes are _monotonic_, i.e. all oriented towards the same quadrant.

This speeds up polygon set operations and boolean ops tremendously, since it allows us to skip a lot of the expensive `O(n^2)` operations.

## Example
=#

