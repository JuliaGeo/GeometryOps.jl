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
# - Planar crossings: validated arithmetic propagates an absolute error bound
#   through endpoint differences, the two determinants, division and coordinate
#   recombination. A coordinate is accepted only when its full interval rounds
#   to one `Float64`; an unsupported range or inconclusive certificate uses the
#   exact `Rational{BigInt}` crossing point instead.
# - Spherical crossings: the **certified correctly rounded** unit direction. The
#   emitted `xyz` is `(RN(x₁), RN(x₂), RN(x₃))` for `x = d/‖d‖` the exact unit
#   crossing direction, `RN` round-to-nearest-even — the planar contract in the
#   chart a spherical node already lives in. Fast path: double-double arithmetic
#   where every intermediate carries a rigorous absolute error bound, so the
#   cancellation in each plane normal (`|a0×a1| = sin θ` with ~eps absolute
#   error, a relative error of eps/θ that no gate on the angle BETWEEN the
#   planes can see) and the near-tangency loss in `na×nb` land in the bound
#   instead of in the output; a coordinate is accepted iff its interval provably
#   rounds to one Float64. Fallback: the exact `_sph_crossing_dir`, each
#   component correctly rounded by an exact midpoint test. WHICH of the two
#   antipodal candidates is meant is a decision, settled in dd when the bound
#   allows and by the exact orients otherwise. The lon/lat row applies its trig
#   to this point, uncertified by design — no decision ever consumes an emitted
#   coordinate.
#
#   The contract is what `_ring_is_subgrid` (maximal_edge_ring.jl) is derived
#   from: emission displaces a vertex by at most ½ ulp per coordinate, so
#   `‖emitted − x‖ ≤ 2⁻⁵³`, and `rk_normalize_usp` is the identity on the result
#   (`|‖u‖² − 1| ≤ 2⁻⁵² + 2 eps`, inside its window), so a cascade level's
#   output re-enters the next level bit-for-bit.
#
#   A crossing has no exact Float64 image in EITHER chart — its position is a
#   `Rational{BigInt}` direction with no finite decimal form — so "no rounding"
#   is not on offer here and nothing below claims it. What the `xyz` row buys is
#   rounding ONCE, in the chart the direction already lives in, instead of twice
#   (normalize to Float64 xyz, then trigonometry to lon/lat). The second
#   rounding is the same one the vertex row pays, with the same latitude
#   dependence.

const _CrossingFloats = ValidatedFloats.CrossingFloats
const _CrossingFloat = _CrossingFloats.CrossingFloat

# Fast certified planar crossing of (a0,a1) × (b0,b1). Returns (x, y, certified).
function _certified_crossing(a0, a1, b0, b1)
    ax0, ay0 = Float64(GI.x(a0)), Float64(GI.y(a0))
    ax1, ay1 = Float64(GI.x(a1)), Float64(GI.y(a1))
    bx0, by0 = Float64(GI.x(b0)), Float64(GI.y(b0))
    bx1, by1 = Float64(GI.x(b1)), Float64(GI.y(b1))
    ax0, ay0 = _CrossingFloat(ax0), _CrossingFloat(ay0)
    ax1, ay1 = _CrossingFloat(ax1), _CrossingFloat(ay1)
    bx0, by0 = _CrossingFloat(bx0), _CrossingFloat(by0)
    bx1, by1 = _CrossingFloat(bx1), _CrossingFloat(by1)
    da_x, da_y = ax1 - ax0, ay1 - ay0
    db_x, db_y = bx1 - bx0, by1 - by0
    c0_x, c0_y = bx0 - ax0, by0 - ay0
    den = da_x * db_y - da_y * db_x
    t = (c0_x * db_y - c0_y * db_x) / den
    x = ax0 + t * da_x
    y = ay0 + t * da_y
    xf = _CrossingFloats.certify(Float64, x)
    yf = _CrossingFloats.certify(Float64, y)
    return (something(xf, _CrossingFloats.center(x)),
            something(yf, _CrossingFloats.center(y)),
            !isnothing(xf) & !isnothing(yf))
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

# ## Validated spherical arithmetic
#
# `CrossingFloat` owns the bounded double-double operations. Its standard
# `cross` and `dot` methods keep the calculation below in vector notation while
# preserving the per-operation error bounds needed after cancellation.

# ## Certified spherical crossing
#
# The arithmetic package admits its proven input range and fails closed outside
# it; degenerate crossing directions likewise route to the exact fallback.

# Fast path: `(certified, x, y, z)` with `(x, y, z)` the correctly rounded unit
# crossing direction when certified.
function _certified_sph_crossing_fast(k::NodeKey; safety = 2.0)
    # The package owns its validated margin; legacy overrides fail closed.
    safety == 2.0 || return (false, 0.0, 0.0, 0.0)
    A0 = _vec3(False(), k.pt); A1 = _vec3(False(), k.a1)
    B0 = _vec3(False(), k.b0); B1 = _vec3(False(), k.b1)
    a0 = UnitSphericalPoint(_CrossingFloat.(A0))
    a1 = UnitSphericalPoint(_CrossingFloat.(A1))
    b0 = UnitSphericalPoint(_CrossingFloat.(B0))
    b1 = UnitSphericalPoint(_CrossingFloat.(B1))
    d = cross(cross(a0, a1), cross(b0, b1))
    #-- which antipode: `d = [b0,b1,a0]·a1 − [b0,b1,a1]·a0` with opposite-sign
    #-- orients at a proper crossing, so `sign(d·(a0+a1)) = sign([b0,b1,a0])` —
    #-- the sign `_crossing_dir_is_positive` computes from exact orients, which
    #-- for short arcs fall under ExactPredicates' filter and cost a BigInt
    #-- escalation each. Decided here with a validated dot product when its sign
    #-- is certified.
    tsign = _CrossingFloats.certified_sign(dot(d, a0 + a1))
    pos = isnothing(tsign) ? _crossing_dir_is_positive(k) : tsign > 0
    pos || (d = -d)
    m = max(abs(_CrossingFloats.center(d[1])), abs(_CrossingFloats.center(d[2])),
            abs(_CrossingFloats.center(d[3])))
    (m > 0 && isfinite(m)) || return (false, 0.0, 0.0, 0.0)
    #-- scale so max |hi| ∈ [1, 2): a power of two, so the scaling is exact
    p = _CrossingFloats.PowerOfTwo(-exponent(m))
    d = map(x -> _CrossingFloats.scale_pow2(x, p), d)
    σ = sqrt(dot(d, d))
    x = _CrossingFloats.certify(Float64, d[1] / σ)
    x === nothing && return (false, 0.0, 0.0, 0.0)
    y = _CrossingFloats.certify(Float64, d[2] / σ)
    y === nothing && return (false, 0.0, 0.0, 0.0)
    z = _CrossingFloats.certify(Float64, d[3] / σ)
    z === nothing && return (false, 0.0, 0.0, 0.0)
    return (true, x, y, z)
end

# Exact fallback: the rational direction, each component `dᵢ/√(d·d)` correctly
# rounded by `_round_div_sqrt`.
function _certified_sph_crossing_exact(k::NodeKey)
    d = _sph_crossing_dir(True(), k)
    S = d[1] * d[1] + d[2] * d[2] + d[3] * d[3]
    return (_round_div_sqrt(d[1], S), _round_div_sqrt(d[2], S), _round_div_sqrt(d[3], S))
end

# The emitted point of a spherical crossing node: the componentwise correctly
# rounded exact unit direction, whichever path produced it. `+ 0.0` clears a
# negative zero so both paths agree bit-for-bit.
function _certified_sph_crossing(k::NodeKey)
    ok, x, y, z = _certified_sph_crossing_fast(k)
    ok || ((x, y, z) = _certified_sph_crossing_exact(k))
    return UnitSphericalPoint(x + 0.0, y + 0.0, z + 0.0)
end

# Spherical → unit-sphere xyz (the default). A vertex node's coordinate IS its
# kernel point: it was normalized once at ingest and is emitted unchanged, so an
# uncut input vertex survives an overlay bit-for-bit. A crossing is the exact
# direction correctly rounded to unit length — one rounding, in the chart the
# direction is already expressed in.
function _emit_node_coord(k::NodeKey{P}, ::Type{P}) where {P <: UnitSphericalPoint}
    k.is_crossing || return k.pt
    return _certified_sph_crossing(k)
end

# Spherical → (lon, lat) degrees. Both arms go through the same trigonometry,
# and for a vertex node that is a round trip out of and back into a chart the
# coordinate did not need to leave — the one place in this file where a value
# that has an exact image in the output format does not get it.
function _emit_node_coord(k::NodeKey{<:UnitSphericalPoint}, ::Type{Tuple{Float64, Float64}})
    k.is_crossing || return _usp_to_lonlat(k.pt)
    return _usp_to_lonlat(_certified_sph_crossing(k))
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
