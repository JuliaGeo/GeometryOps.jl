# # Spherical RelateKernel
#
#=
Spherical implementation of the RelateKernel contract declared in `kernel.jl`,
over `UnitSphericalPoint{Float64}`. Every predicate is a sign of
det(u, v, w) = (u×v)·w, so the exact path mirrors planar: a float filter then an
exact fallback (`ExactPredicates.orient` for the plain orient, `Rational{BigInt}`
on the xyz components for composites). No intersection coordinate is ever
constructed. See the design doc 2026-06-15.

All of `cross`, `⋅`, `normalize` (LinearAlgebra), `ExactPredicates`, and the
`UnitSpherical` names are already in scope here — this file is `include`d into
`GeometryOps`, which `using`s them at the top of the module. The
`_rebuild_point(::UnitSphericalPoint, …)` hook that keeps node points typed
lives next to the generic `_rebuild_point` in `kernel.jl`.
=#

# xyz tuple of a 3D point, for the ExactPredicates / Rational{BigInt} paths.
@inline _tup3(u) = (GI.x(u), GI.y(u), GI.z(u))
