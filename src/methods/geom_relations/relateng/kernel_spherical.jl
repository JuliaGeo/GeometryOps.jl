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

# ## rk_orient

# Orientation of `c` relative to the great-circle arc `(a, b)`: the sign of the
# scalar triple product (a×b)·c. Exact path: `ExactPredicates.orient` over the
# xyz tuples about the origin (the spike measured 3.5 ns); NOT
# `UnitSpherical.spherical_orient`, whose eps*16 tolerance is unfit for the
# exact contract. Float path: the plain triple product.
rk_orient(::Spherical, a, b, c; exact) = _rk_orient(booltype(exact), a, b, c)
@inline _rk_orient(::True, a, b, c) =
    ExactPredicates.orient(_tup3(a), _tup3(b), _tup3(c), (0.0, 0.0, 0.0))
@inline _rk_orient(::False, a, b, c) = cross(a, b) ⋅ c

# ## rk_point_on_segment

# Whether `p` lies on the closed minor arc `[q0, q1]`. Two conditions: `p` is on
# the arc's great circle (coplanar with `q0`, `q1`, the origin — an exact orient
# == 0), and within the minor-arc span. The span test uses dots as cos-of-angle:
# for unit vectors `q0·p >= q0·q1` means `p` is no farther from `q0` than `q1`
# is, and symmetrically for `q1`; together they pin `p` to the minor arc. Engine
# inputs are unit; the conformance grid points are chosen so the dot signs are
# exact.
function rk_point_on_segment(m::Spherical, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    qq = q0 ⋅ q1
    return (q0 ⋅ p) >= qq && (q1 ⋅ p) >= qq
end
