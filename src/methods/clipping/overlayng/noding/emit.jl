# # Emission (design §2.6) — the only lossy step
#
# `node_point(arr, id)` realizes a node's OUTPUT coordinate on demand and caches
# it. Nothing in the arrangement's decisions ever consumed a constructed
# coordinate; this is where the exact symbolic result is rounded to Float64 for
# output. Both manifolds emit `(x, y)` / `(lon, lat)` as `Tuple{Float64,Float64}`.
#
# - Vertex nodes: the input vertex, bit-exact pass-through.
# - Planar crossings: a **certified double-double** fast path (spike S3, 100%
#   certified on 64,982 real crossings, 0 disagreements with the rational answer,
#   273×) — TwoSum on endpoint differences, compensated 2×2 determinants,
#   dd division, dd recombination; the coordinate is accepted iff its residual
#   plus the dd error bound is below ½ ulp, with the determinant-conditioning
#   term so near-parallel pairs fail the certificate. Fallback: the exact
#   `Rational{BigInt}` crossing point, rounded.
# - Spherical crossings: the Float64 crossing direction `±(na×nb)`, accepted when
#   the arcs clear a near-tangency conditioning gate (spike S3 measured the float
#   direction at ≤1.4e-14° ≈ 1.5 nm; the lon/lat trig itself is uncertified by
#   design — no decision ever consumes an emitted coordinate). Fallback: the
#   exact `_sph_crossing_dir`, normalized and converted.

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

# ## Per-manifold node coordinate realization (dispatched on the kernel point `P`)

# Planar: vertex pass-through; crossing via the certified dd path with a rational
# fallback (identical to the fallback the S3 audit compared against).
function _emit_node_coord(k::NodeKey{Tuple{Float64, Float64}})
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

function _emit_node_coord(k::NodeKey{<:UnitSphericalPoint})
    k.is_crossing || return _usp_to_lonlat(k.pt)
    #-- float na, nb, d = na × nb, with the conditioning gate
    A0 = _vec3(False(), k.pt); A1 = _vec3(False(), k.a1)
    B0 = _vec3(False(), k.b0); B1 = _vec3(False(), k.b1)
    na = _cross3(A0, A1); nb = _cross3(B0, B1)
    d = _cross3(na, nb)
    d2 = _dot3(d, d); na2 = _dot3(na, na); nb2 = _dot3(nb, nb)
    if d2 >= _SPH_TANGENT_GATE^2 * na2 * nb2
        dir = _sph_crossing_dir(False(), k)          # picks the interior candidate
        return _dir_to_lonlat(dir)
    end
    #-- near-tangent fallback: exact direction (Rational), then normalize + convert
    return _dir_to_lonlat(_sph_crossing_dir(True(), k))
end

@inline function _dir_to_lonlat(d)
    s = sqrt(Float64(d[1])^2 + Float64(d[2])^2 + Float64(d[3])^2)
    return _usp_to_lonlat(UnitSphericalPoint(Float64(d[1]) / s, Float64(d[2]) / s, Float64(d[3]) / s))
end

@inline function _usp_to_lonlat(u)
    ll = GeographicFromUnitSphere()(u)
    return (Float64(ll[1]), Float64(ll[2]))
end

# ## Public-internal accessor

"""
    node_point(arr::NodedArrangement, id) -> Tuple{Float64,Float64}

The realized output coordinate of node `id` (planar `(x, y)` / spherical
`(lon, lat)`), memoized in the node table (design §2.6). The only place a
constructed coordinate enters the substrate.
"""
function node_point(arr::NodedArrangement, id::Integer)
    t = arr.nodes
    i = Int(id)
    @inbounds t.realized[i] && return t.coords[i]
    c = _emit_node_coord(t.keys[i])
    @inbounds t.coords[i] = c
    @inbounds t.realized[i] = true
    return c
end

# The two output coordinates of a noded edge (convenience for callers/tests).
edge_endpoints(arr::NodedArrangement, e::NodedEdge) =
    (node_point(arr, e.node_lo), node_point(arr, e.node_hi))
