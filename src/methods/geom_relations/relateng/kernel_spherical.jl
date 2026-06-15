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

# ## Exact-aware 3-vector arithmetic
#
# Composite predicates (arc membership, proper crossing, node coincidence)
# reduce to signs of polynomials in the xyz components. With `exact = True()` we
# evaluate over `Rational{BigInt}` (Float64 are dyadic rationals → exact); with
# `False()`, Float64. `_vec3(bt, p)` lifts a point to the chosen number type; the
# rest are plain tuple cross/dot, so one code path serves both — exactly how the
# planar kernel threads `exact`.
@inline _vec3(::True, u) = (Rational{BigInt}(GI.x(u)), Rational{BigInt}(GI.y(u)), Rational{BigInt}(GI.z(u)))
@inline _vec3(::False, u) = (Float64(GI.x(u)), Float64(GI.y(u)), Float64(GI.z(u)))
@inline _cross3(a, b) = (a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
@inline _dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline _iszero3(a) = iszero(a[1]) && iszero(a[2]) && iszero(a[3])
@inline _neg3(a) = (-a[1], -a[2], -a[3])
# `w` strictly interior to the minor arc (a, b) with normal n = a×b.
@inline _strictly_in_arc3(w, a, b, n) = _dot3(_cross3(a, w), n) > 0 && _dot3(_cross3(w, b), n) > 0
_usp_eq(p, q) = GI.x(p) == GI.x(q) && GI.y(p) == GI.y(q) && GI.z(p) == GI.z(q)

# ## rk_point_on_segment

# Whether `p` lies on the closed minor arc `[q0, q1]`. Two conditions: `p` is on
# the arc's great circle (coplanar with `q0`, `q1`, origin — an exact orient ==
# 0), and within the minor-arc span. The span test is scale-invariant: writing
# the coplanar `p` as `α q0 + β q1`, `p` is on the closed minor arc iff α, β ≥ 0,
# and `sign(β) = sign((q0×p)·n)`, `sign(α) = sign((p×q1)·n)` with `n = q0×q1` — a
# pure determinant sign, correct for unit and non-unit inputs alike.
function rk_point_on_segment(m::Spherical, p, q0, q1; exact)
    rk_orient(m, q0, q1, p; exact) == 0 || return false
    return _on_arc_span(booltype(exact), p, q0, q1)
end
@inline function _on_arc_span(bt, p, q0, q1)
    P = _vec3(bt, p); Q0 = _vec3(bt, q0); Q1 = _vec3(bt, q1)
    n = _cross3(Q0, Q1)
    return _dot3(_cross3(Q0, P), n) >= 0 && _dot3(_cross3(P, Q1), n) >= 0
end

# ## Ingest and interaction bounds

# Renormalize to unit length (Float32-sourced data — e.g. Natural Earth GeoJSON
# converted to Float64 — is ~1e-8 off unit and trips `robust_cross_product`).
@inline rk_normalize_usp(u) = UnitSphericalPoint(normalize(u))

# Canonical kernel point of a GeoInterface point: lon/lat (2D) → unit xyz, or an
# already-3D point treated as xyz; renormalized and signed-zero normalized so
# the same vertex always produces identical bits (NodeKey equality). The vertex
# ingest (Phase 3 `_to_kernel_point`) and the extent computation below share
# this, so a vertex and its extent agree exactly.
@inline function _spherical_kernel_point(p)
    u = GI.is3d(p) ?
        UnitSphericalPoint(Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
        UnitSphereFromGeographic()((Float64(GI.x(p)), Float64(GI.y(p))))
    return _node_point(rk_normalize_usp(u))
end

# 3D AABB of a minor great-circle arc (spike-proven, 0/102k fuzz escapes). A
# great-circle arc bulges outside the coordinate box of its endpoints; the
# extremum of coordinate i on the circle with normal n is at ±w, the normalized
# projection of axis eᵢ onto the circle's plane. The box is extended by whichever
# of ±w lies on the minor arc; a few ulps of widening absorb the roundoff in w.
@inline _on_minor_arc(w, a, b, n) = (cross(a, w) ⋅ n) >= 0.0 && (cross(w, b) ⋅ n) >= 0.0
@inline _widen(lo, hi) = (prevfloat(lo, 4), nextfloat(hi, 4))

function arc_extent(a, b)
    n = cross(a, b)
    n2 = n ⋅ n
    xlo, xhi = minmax(a[1], b[1]); ylo, yhi = minmax(a[2], b[2]); zlo, zhi = minmax(a[3], b[3])
    if n2 > 0.0
        invn2 = inv(n2)
        @inbounds for i in 1:3
            ei_n = n[i]
            wx = (i == 1) - ei_n * n[1] * invn2
            wy = (i == 2) - ei_n * n[2] * invn2
            wz = (i == 3) - ei_n * n[3] * invn2
            wnorm = sqrt(wx^2 + wy^2 + wz^2)
            wnorm == 0.0 && continue   # axis ⟂ plane: coordinate constant 0, endpoints cover it
            w = UnitSphericalPoint(wx / wnorm, wy / wnorm, wz / wnorm)
            for ww in (w, -w)
                if _on_minor_arc(ww, a, b, n)
                    ci = ww[i]
                    if i == 1; xlo = min(xlo, ci); xhi = max(xhi, ci)
                    elseif i == 2; ylo = min(ylo, ci); yhi = max(yhi, ci)
                    else; zlo = min(zlo, ci); zhi = max(zhi, ci)
                    end
                end
            end
        end
    end
    return Extents.Extent(X = _widen(xlo, xhi), Y = _widen(ylo, yhi), Z = _widen(zlo, zhi))
end

@inline _point_box(u) =
    Extents.Extent(X = _widen(GI.x(u), GI.x(u)), Y = _widen(GI.y(u), GI.y(u)), Z = _widen(GI.z(u), GI.z(u)))

# ## rk_classify_intersection
#
# The two great circles meet at ±d, d = (a0×a1)×(b0×b1). `SS_PROPER` iff one of
# ±d is strictly interior to both minor arcs (the candidate-direct formulation —
# planar straddle tests are NOT sufficient on the sphere, where arcs can straddle
# each other's great circle while meeting only at the antipodal point). Endpoint
# incidences are exact arc-membership; collinear = the arcs share a great circle
# (d == 0). No intersection coordinate is constructed.
function rk_classify_intersection(m::Spherical, a0, a1, b0, b1; exact)
    a0_on_b = rk_point_on_segment(m, a0, b0, b1; exact)
    a1_on_b = rk_point_on_segment(m, a1, b0, b1; exact)
    b0_on_a = rk_point_on_segment(m, b0, a0, a1; exact)
    b1_on_a = rk_point_on_segment(m, b1, a0, a1; exact)
    return _sph_classify(booltype(exact), a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
end

function _sph_classify(bt, a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    A0 = _vec3(bt, a0); A1 = _vec3(bt, a1); B0 = _vec3(bt, b0); B1 = _vec3(bt, b1)
    na = _cross3(A0, A1); nb = _cross3(B0, B1)
    d = _cross3(na, nb)
    n_inc = a0_on_b + a1_on_b + b0_on_a + b1_on_a
    if _iszero3(d)   # same great circle (or a degenerate, zero-length arc)
        n_inc == 0 && return SegSegClass(SS_DISJOINT, false, false, false, false)
        # a degenerate (zero-length) arc on the other is a touch, not an overlap
        zero_len = _iszero3(na) || _iszero3(nb)
        shared_only = n_inc == 2 && (a0_on_b || a1_on_b) && (b0_on_a || b1_on_a) &&
            (_usp_eq(a0, b0) || _usp_eq(a0, b1) || _usp_eq(a1, b0) || _usp_eq(a1, b1))
        kind = (shared_only || zero_len) ? SS_TOUCH : SS_COLLINEAR
        return SegSegClass(kind, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    if a0_on_b || a1_on_b || b0_on_a || b1_on_a
        return SegSegClass(SS_TOUCH, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    end
    nd = _neg3(d)
    if (_strictly_in_arc3(d, A0, A1, na) && _strictly_in_arc3(d, B0, B1, nb)) ||
       (_strictly_in_arc3(nd, A0, A1, na) && _strictly_in_arc3(nd, B0, B1, nb))
        return SegSegClass(SS_PROPER, false, false, false, false)
    end
    return SegSegClass(SS_DISJOINT, false, false, false, false)
end

# Interaction bounds on the sphere: a 3D `Extent{(:X,:Y,:Z)}` in unit-sphere xyz
# (the engine works in xyz after ingest), as the union of `arc_extent` over the
# geometry's edges. Area-element interiors reach beyond their boundary slab — the
# ±eᵢ axis-point extension is added in Task 11.
rk_interaction_bounds(m::Spherical, geom) = _sph_bounds(m, GI.trait(geom), geom)

function _sph_bounds(::Spherical, ::GI.AbstractPointTrait, geom)
    return _point_box(_spherical_kernel_point(geom))
end
function _sph_bounds(::Spherical, ::GI.AbstractCurveTrait, geom)
    n = GI.npoint(geom)
    prev = _spherical_kernel_point(GI.getpoint(geom, 1))
    n == 1 && return _point_box(prev)
    ext = nothing
    for i in 2:n
        cur = _spherical_kernel_point(GI.getpoint(geom, i))
        e = arc_extent(prev, cur)
        ext = ext === nothing ? e : Extents.union(ext, e)
        prev = cur
    end
    return ext
end
function _sph_bounds(m::Spherical, ::GI.AbstractPolygonTrait, geom)
    ext = _sph_bounds(m, GI.trait(GI.getexterior(geom)), GI.getexterior(geom))
    for hole in GI.gethole(geom)
        ext = Extents.union(ext, _sph_bounds(m, GI.trait(hole), hole))
    end
    return ext
end
function _sph_bounds(m::Spherical, ::GI.AbstractGeometryTrait, geom)
    ext = nothing
    for g in GI.getgeom(geom)
        e = _sph_bounds(m, GI.trait(g), g)
        ext = ext === nothing ? e : Extents.union(ext, e)
    end
    return ext
end
