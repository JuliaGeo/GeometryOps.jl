# NOTE: This functionality is experimental and may change at any time.

# # Emission (design §2.6) — the only lossy step
#
# `node_point(arr, id)` realizes a node's OUTPUT coordinate on demand and caches
# it. Nothing in the arrangement's decisions ever consumed a constructed
# coordinate; this is where the exact symbolic result is rounded to Float64 for
# output.
#
# The output point type is a parameter of the arrangement, not a property of the
# manifold (`OverlayNG`'s `point_type`), so this file has one method per
# (kernel point, output point) pair. There are three:
#
#   planar    → `Tuple{Float64,Float64}`         `(x, y)`
#   spherical → `UnitSphericalPoint{Float64}`    unit-sphere `xyz` (the default)
#   spherical → `Tuple{Float64,Float64}`         `(lon, lat)` in degrees
#
# - Vertex nodes: the input vertex, bit-exact pass-through — on the two rows
#   whose output type IS the kernel type this is literally `k.pt`. The spherical
#   lon/lat row is the exception, and the only one: it sends a coordinate that
#   already has an exact image in the output format back out through
#   `atan`/`asin`. Measured over 200 000 uniformly random directions that round
#   trip displaces the point by ≤3.2 ULPs of the unit sphere below 60° latitude,
#   7.9 at 75° and 126 above 89.5°; along single parallels, 1 364 ULPs at 89.99°
#   and 2 915 (6.5e-13 rad) at 89.999°. The blow-up is the chart's — a degree of
#   longitude is `cos φ` of an arc — and it is what the default row removes.
# - Planar crossings: a **certified double-double** fast path (spike S3, 100%
#   certified on 64,982 real crossings, 0 disagreements with the rational answer,
#   273×) — TwoSum on endpoint differences, compensated 2×2 determinants,
#   dd division, dd recombination; the coordinate is accepted iff its residual
#   plus the dd error bound is below ½ ulp, with the determinant-conditioning
#   term so near-parallel pairs fail the certificate. Fallback: the exact
#   `Rational{BigInt}` crossing point, rounded.
# - Spherical crossings: the Float64 crossing direction `±(na×nb)`, accepted when
#   the arcs clear a near-tangency conditioning gate (spike S3 measured the float
#   direction at ≤1.4e-14° ≈ 1.5 nm; the trig on the lon/lat row is uncertified
#   by design — no decision ever consumes an emitted coordinate). Fallback: the
#   exact `_sph_crossing_dir`, normalized and converted. The gate bounds the
#   direction's *magnitude* error only; WHICH of the two antipodal candidates is
#   meant is a decision and `_sph_crossing_dir` settles it exactly on both paths.
#
#   A crossing has no exact Float64 image in EITHER chart — its position is a
#   `Rational{BigInt}` direction with no finite decimal form — so "no rounding"
#   is not on offer here and nothing below claims it. What the `xyz` row buys is
#   rounding ONCE, in the chart the direction already lives in, instead of twice
#   (normalize to Float64 xyz, then trigonometry to lon/lat). The second
#   rounding is the same one the vertex row pays, with the same latitude
#   dependence.

# ## Error-free transforms and double-double primitives (spike S3, productionized)

@inline function _twosum(a::Float64, b::Float64)
    s = a + b; bb = s - a
    return (s, (a - (s - bb)) + (b - bb))
end
@inline function _twoproduct(a::Float64, b::Float64)
    p = a * b
    return (p, fma(a, b, -p))
end
@inline _diff_dd(a::Float64, b::Float64) = _twosum(a, -b)   # (hi, lo) == a - b exactly

@inline function _ddmul(ah, al, bh, bl)                     # (ah+al)*(bh+bl)
    (ph, pl) = _twoproduct(ah, bh)
    return _twosum(ph, pl + (ah * bl + al * bh))
end
@inline function _ddsub(ah, al, bh, bl)                     # (ah+al) - (bh+bl)
    (sh, se) = _twosum(ah, -bh)
    return _twosum(sh, (se + al) - bl)
end
@inline function _ddadd(ah, al, bh, bl)                     # (ah+al) + (bh+bl)
    (sh, se) = _twosum(ah, bh)
    return _twosum(sh, (se + al) + bl)
end
# det = a*d - b*c, each operand a double-double 2-tuple
@inline function _det2_ddfull(a, b, c, d)
    (mh, ml) = _ddmul(a[1], a[2], d[1], d[2])
    (nh, nl) = _ddmul(b[1], b[2], c[1], c[2])
    return _ddsub(mh, ml, nh, nl)
end
# dd / dd -> dd (Dekker)
@inline function _div_dd(ah, al, bh, bl)
    q1 = ah / bh
    (ph, pl) = _twoproduct(q1, bh)
    (sh, sl) = _twosum(ah, -ph)
    r = ((sh - pl) + sl) + al - q1 * bl
    return _twosum(q1, r / bh)
end

# Certified correctly-rounded emit of one coordinate: `xf = fl(hi+lo)`, exact
# residual `rem` from TwoSum, accept iff `|rem| + dderr < ½ ulp(xf)`.
@inline function _certify_coord(hi, lo, dderr)
    (xf, rem) = _twosum(hi, lo)
    return (xf, abs(rem) + dderr < 0.5 * eps(xf))
end

# Fast certified planar crossing of (a0,a1) × (b0,b1). Returns (x, y, certified).
function _certified_crossing(a0, a1, b0, b1)
    ax0, ay0 = Float64(GI.x(a0)), Float64(GI.y(a0))
    ax1, ay1 = Float64(GI.x(a1)), Float64(GI.y(a1))
    bx0, by0 = Float64(GI.x(b0)), Float64(GI.y(b0))
    bx1, by1 = Float64(GI.x(b1)), Float64(GI.y(b1))
    #-- exact endpoint differences as double-doubles
    da_x = _diff_dd(ax1, ax0); da_y = _diff_dd(ay1, ay0)
    db_x = _diff_dd(bx1, bx0); db_y = _diff_dd(by1, by0)
    c0_x = _diff_dd(bx0, ax0); c0_y = _diff_dd(by0, ay0)
    (denh, denl) = _det2_ddfull(da_x, da_y, db_x, db_y)   # da × db
    (tnh, tnl)   = _det2_ddfull(c0_x, c0_y, db_x, db_y)   # (b0-a0) × db
    (th, tl) = _div_dd(tnh, tnl, denh, denl)
    (txh, txl) = _ddmul(th, tl, da_x[1], da_x[2]); (xh, xl) = _ddadd(ax0, 0.0, txh, txl)
    (tyh, tyl) = _ddmul(th, tl, da_y[1], da_y[2]); (yh, yl) = _ddadd(ay0, 0.0, tyh, tyl)
    #-- dd error bounds amplified by determinant conditioning: near-parallel ⇒
    #-- small |denom| ⇒ large condK ⇒ certificate correctly fails (spike S3)
    u2 = eps(Float64)^2
    scale = abs(da_x[1]) + abs(da_y[1]) + abs(db_x[1]) + abs(db_y[1])
    condK = (abs(da_x[1] * db_y[1]) + abs(da_y[1] * db_x[1])) / max(abs(denh), floatmin(Float64))
    tmag = abs(th)
    ex = 64 * u2 * (abs(xh) + tmag * abs(da_x[1]) * condK + scale)
    ey = 64 * u2 * (abs(yh) + tmag * abs(da_y[1]) * condK + scale)
    (xf, cx) = _certify_coord(xh, xl, ex)
    (yf, cy) = _certify_coord(yh, yl, ey)
    return (xf, yf, cx & cy)
end

# ## Node coordinate realization (dispatched on the kernel point AND the output type)

# Planar: vertex pass-through; crossing via the certified dd path with a rational
# fallback (identical to the fallback the S3 audit compared against).
function _emit_node_coord(k::NodeKey{Tuple{Float64, Float64}},
        ::Type{Tuple{Float64, Float64}})
    k.is_crossing || return k.pt
    (x, y, cert) = _certified_crossing(k.pt, k.a1, k.b0, k.b1)
    cert && return (x, y)
    rx, ry = _exact_crossing_point(k)
    return (Float64(rx), Float64(ry))
end

# Near-tangency gate for the spherical float direction: |na×nb|² ≥ tol²·|na|²·|nb|²
# means the arcs' planes meet at ≥ ~1e-9 rad, so the float direction's relative
# error (≈ eps / sin θ) is bounded well below the ≤1.4e-14° the design accepts.
# Below the gate the crossing is near-tangent and falls to the exact direction.
const _SPH_TANGENT_GATE = 1e-9

# The crossing direction the two spherical rows share: the gated Float64
# `±(na × nb)` when the arcs' planes are far enough from parallel, the exact
# rational direction otherwise. Unnormalized on both branches — each row
# normalizes into its own output chart, so nothing is normalized twice.
function _sph_emit_dir(k::NodeKey{<:UnitSphericalPoint})
    #-- float na, nb, d = na × nb, with the conditioning gate
    A0 = _vec3(False(), k.pt); A1 = _vec3(False(), k.a1)
    B0 = _vec3(False(), k.b0); B1 = _vec3(False(), k.b1)
    na = _cross3(A0, A1); nb = _cross3(B0, B1)
    d = _cross3(na, nb)
    d2 = _dot3(d, d); na2 = _dot3(na, na); nb2 = _dot3(nb, nb)
    #-- `_sph_crossing_dir` picks the interior candidate of the antipodal pair
    d2 >= _SPH_TANGENT_GATE^2 * na2 * nb2 && return _sph_crossing_dir(False(), k)
    #-- near-tangent fallback: the exact direction (Rational)
    return _sph_crossing_dir(True(), k)
end

# Spherical → unit-sphere xyz (the default). A vertex node's coordinate IS its
# kernel point: it was normalized once at ingest and is emitted unchanged, so an
# uncut input vertex survives an overlay bit-for-bit. A crossing is the exact
# direction normalized to unit length — one rounding, in the chart the direction
# is already expressed in.
function _emit_node_coord(k::NodeKey{P}, ::Type{P}) where {P <: UnitSphericalPoint}
    k.is_crossing || return k.pt
    return _dir_to_usp(_sph_emit_dir(k))
end

# Spherical → (lon, lat) degrees. Both arms go through the same trigonometry,
# and for a vertex node that is a round trip out of and back into a chart the
# coordinate did not need to leave — the one place in this file where a value
# that has an exact image in the output format does not get it.
function _emit_node_coord(k::NodeKey{<:UnitSphericalPoint}, ::Type{Tuple{Float64, Float64}})
    k.is_crossing || return _usp_to_lonlat(k.pt)
    return _dir_to_lonlat(_sph_emit_dir(k))
end

# A crossing direction (Float64 or `Rational{BigInt}` components) as a unit
# `UnitSphericalPoint{Float64}` — the same construction `_node_kernel_point`
# uses for the exact kernel position of a crossing node, so the emitted point
# and the kernel point of a crossing are built by one formula.
@inline _dir_to_usp(d) =
    rk_normalize_usp(UnitSphericalPoint(Float64(d[1]), Float64(d[2]), Float64(d[3])))

@inline function _dir_to_lonlat(d)
    s = sqrt(Float64(d[1])^2 + Float64(d[2])^2 + Float64(d[3])^2)
    return _usp_to_lonlat(UnitSphericalPoint(Float64(d[1]) / s, Float64(d[2]) / s, Float64(d[3]) / s))
end

@inline function _usp_to_lonlat(u)
    ll = GeographicFromUnitSphere()(u)
    return (Float64(ll[1]), Float64(ll[2]))
end

# ## Properties of the output point type
#
# The engine asks the output format three things, and they are all answered
# here, beside the emitter, so that adding a fourth output type is one block of
# edits rather than a hunt: how many numbers a coordinate has (ring bounding
# boxes), whether GeoInterface should call the result 3D, and — in
# maximal_edge_ring.jl, which needs the ring's own location to answer it — how
# far apart the format's representable positions are.

@inline _output_ncoord(::Type{Tuple{Float64, Float64}}) = 2
@inline _output_ncoord(::Type{<:UnitSphericalPoint}) = 3

# One (min, max) pair per coordinate.
@inline _bbox_len(::Type{T}) where {T} = 2 * _output_ncoord(T)

@inline _output_is3d(::Type{Tuple{Float64, Float64}}) = false
@inline _output_is3d(::Type{<:UnitSphericalPoint}) = true

# ## Public-internal accessor

"""
    node_point(arr::NodedArrangement{P,T}, id) -> T

The realized output coordinate of node `id`, in the arrangement's output point
type `T` (see `output_point_type`), memoized in the node table (design §2.6).
The only place a constructed coordinate enters the substrate.
"""
function node_point(arr::NodedArrangement{P, T}, id::Integer) where {P, T}
    t = arr.nodes
    i = Int(id)
    @inbounds t.realized[i] && return t.coords[i]
    c = _emit_node_coord(t.keys[i], T)::T
    @inbounds t.coords[i] = c
    @inbounds t.realized[i] = true
    return c
end

# The two output coordinates of a noded edge (convenience for callers/tests).
edge_endpoints(arr::NodedArrangement, e::NodedEdge) =
    (node_point(arr, e.node_lo), node_point(arr, e.node_hi))
