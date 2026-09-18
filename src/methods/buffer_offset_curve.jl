# # Raw offset curves

#=
The linework half of the native [`buffer`](@ref): given a geometry and a signed
distance, produce the closed rings whose winding-number-`>= 1` region is the
buffer. [`_winding_overlay`](@ref) does the region half.

## What a raw offset curve is

Walk a ring in its canonical direction and, at every vertex, step across to the
curve offset by `|distance|`:

- **outside turn** — the offset points on either side of the vertex are joined by
  a circular fillet around it, approximated by `quadsegs` segments per quadrant,
  or by a single point when they are closer together than `0.05·|d|`;
- **inside turn** — the two offset points are joined *through the vertex itself*,
  which traces a small loop, unless the guarded trim below replaces the loop by
  the single point where the two offset segments cross;
- **collinear** — nothing is emitted (the two offsets are collinear too), except
  at an exact reversal, which gets a 180° fillet around the tip.

The result is one closed, freely self-intersecting ring per input ring, line or
point. The inside-turn loops are traversed in the opposite sense to the rest of
the curve, so they subtract themselves from the winding number and need no
explicit trimming — which is why this generator can be a straight transcription
of JTS `OffsetSegmentGenerator` *minus* its heuristics.

## The guarded inside-turn trim

Every inside-turn loop is one more self-crossing for the noder to find, and the
noder is where the whole operation's time goes. JTS drops the loop (`addInsideTurn`,
`OffsetSegmentGenerator.java:314`) by emitting only the point `X` where the two
offset *segments* meet. That is not free — it deletes the kite `(X, o0.p1, s1,
o1.p0)` from the curve's 1-chain, which shifts the winding number by one inside
the kite — so the trim is taken only under three guards:

1. **The crossing must be proper**: `X` strictly interior to *both* offset
   segments, by exact orientation predicates. Then the kite splits along the
   diagonal `X–s1` into `T0 = (X, o0.p1, s1) ⊆ HR_0` and `T1 = (X, s1, o1.p0) ⊆
   HR_1`, the two half-rectangles the adjacent edges contribute to the coverage
   count, and `∂K = ∂T0 + ∂T1` exactly. The winding field therefore stays of the
   form `1_P ± cov` with `cov` the coverage of the *same* pieces minus at most
   one unit each on `T0` and `T1`, so it can never become a field of the wrong
   shape: a trim can only shrink a dilation and only grow an erosion, and only
   within `|d|` of the input.
2. **No backtracking along an offset segment**: the trim point must lie strictly
   beyond the point the curve already sits at on that segment. Two sharp inside
   turns joined by an edge shorter than the two trims can otherwise remove
   overlapping triangles from the *same* half-rectangle, which is the one way
   the coverage count can go negative and a dilation lose a real piece.
3. **One untrimmed vertex path per side**: the first inside turn of each curve
   side always takes the vertex path. A winding shift only accumulates where a
   point lies in the kites of *every* corner of a ring — the collapse cases, a
   shell eroded past its inradius or a hole shrunk to nothing — and one
   untrimmed corner breaks that cycle. It also fixes the closure: the ring's
   first emitted point is then the end of its first offset segment, so the
   closing edge cannot run backwards either.

When both edges at a trimmed corner are at least `|d|` long the kept region is
not merely close but *identical*. Write `(a, b)` for `HR_1`'s own frame — `a`
along the outgoing edge, `b` toward the offset side — and `θ` for the turn
angle. `T0`'s three vertices are `(0, 0)`, `(|d| sin θ, |d| cos θ)` and
`(|d| tan(θ/2), |d|)`, so on `T0` the linear coordinate `a` stays within
`[0, L₁]` once `L₁ ≥ |d| ≥ |d| sin θ` and `L₁ > |d| tan(θ/2)` (guard 1 gives the
second), and `b ≤ |d|`. The part of `T0` with `b < 0` is on the far side of the
outgoing edge, where the identity's `1_P` term absorbs the shift; the rest lies
in `HR_1`, so every point of the kite keeps a second covering piece and the
threshold `w ≥ 1` cannot move. `T1` follows by the mirror argument. Shorter
edges have to borrow their second cover from further pieces, which is why the
guards above and the property tests, not this argument, are what the trim rests
on.

## Relationship to JTS

The primitives below (`_offset_segment`, `_offset_directed_fillet!`,
`_offset_corner_fillet!`, `_offset_circle`, `_offset_add_point!`) reproduce
`OffsetSegmentGenerator` and `OffsetSegmentString` bit for bit, including
`Angle.cosSnap`'s `5e-16` snap, `addDirectedFillet`'s `(int)(angle/quantum + 0.5)`
segment count (a *rounded*, not ceiled, count — GEOS 3.14.1 does the same), and
`createCircle`'s clockwise sweep from `(x + d, y)`. Buffers therefore agree with
GEOS vertex for vertex on the fillets.

The winding overlay makes three more JTS heuristics unnecessary, so the generator
omits them:

| JTS heuristic | Why it is dropped |
|:--|:--|
| the `k = 80` closing segment | its `(m0, s1, m1)` loop shifts the winding with no kite to bound it, and measured 14/1692 wrong collapse cases even behind the guards above |
| `BufferInputLineSimplifier` | deletes input vertices up to `0.003·d` off their chord; we keep them |
| `isRingCurveInverted` | only needed because a curve trimmed without the guards above can invert |

The first changes the answer by a bounded amount in GEOS's favour of speed
(measured: `<= 2.4e-8` relative area on a 0.002·d concavity, `~1e-10` on
1000-vertex inputs).

`OFFSET_SEGMENT_SEPARATION_FACTOR` *is* kept, and unlike the trim it is an
approximation rather than an identity. At an outside turn of `θ` the offset
endpoints are `2|d| sin(θ/2)` apart, so at GEOS's `1e-3` the shortcut fires
below about 0.06°, where the fillet emits no arc point anyway and the join is
two points a straight edge can replace. What it removes from the chain is the
thin triangle `(o0.p1, o1.p0, o1.p1)`, which lies in `W_v ∪ HR_1` — so the
winding field stays a coverage field and only ever loses, never gains — and it
displaces the boundary by at most `|o0.p1 − o1.p0| ≤ 1e-3·|d|` near that one
vertex. The gain is one output vertex per smooth vertex instead of two: a
10 000-gon buffered outwards drops from 20 000 raw segments to 10 000, which is
also what makes the output vertex count match GEOS's.

## Input cleaning

Mirrors `BufferCurveSetBuilder`: consecutive duplicate points are dropped, a ring
with fewer than three distinct vertices buffers as a line, a line with one
distinct point buffers as a point, and for `distance < 0` only areal parts
contribute at all. Polygons are canonicalised to shell-CCW / holes-CW, so the
material interior is always on the left of travel and the offset side is a
function of the sign of the distance alone.

`isRingFullyEroded` is kept too, as a speed guard: the winding identity gets a
collapsed ring right without it, and pays Θ(n²) self-crossings to do so. Its
soundness is argued where `_ring_fully_eroded` is defined, below, and
`prune_eroded_rings = false` turns it off so the engine's own answer at collapse
can be read directly.
=#

# `OFFSET_SEGMENT_SEPARATION_FACTOR`: at an outside turn whose two offset
# endpoints are closer than this fraction of the distance, the join is one point.
# GEOS 3.14.1 uses `1e-3` where JTS 1.20 uses `0.05`, and this value is what
# keeps the output's vertex count equal to GEOS's: measured, both implementations
# keep two points at a 1000-gon's 6.3e-3·|d| gap and one at a 10 000-gon's 6.3e-4.
const _OFFSET_SEPARATION_FACTOR = 1e-3

# Offset side, in the sense of JTS `Position`: LEFT of the direction of travel.
const _OFFSET_LEFT = 1
const _OFFSET_RIGHT = -1
# Fillet sweep direction, in the sense of JTS `Orientation`.
const _OFFSET_CW = -1
const _OFFSET_CCW = 1

# JTS `Angle.cosSnap` / `sinSnap`: a trig result this close to an axis came from
# an angle that was meant to be axis-aligned, and the fillet vertex should be too.
@inline _offset_snap(x::Float64) = abs(x) < 5e-16 ? 0.0 : x

#=
The segment `(a, b)` displaced by `side * d`, as its two endpoints (JTS
`computeOffsetSegment`). The `sqrt(dx^2 + dy^2)` spelling is JTS's, not `hypot`:
`hypot` is more accurate, and that accuracy would move vertices off GEOS's.
=#
@inline function _offset_segment(a::P, b::P, side::Int, d::Float64) where {P}
    dx = b[1] - a[1]; dy = b[2] - a[2]
    len = sqrt(dx * dx + dy * dy)
    ux = side * d * dx / len; uy = side * d * dy / len
    return (a[1] - uy, a[2] + ux), (b[1] - uy, b[2] + ux)
end

#=
The point list a curve is accumulated into (JTS `OffsetSegmentString`). Works in
the floating precision model, like the winding overlay it feeds, and keeps JTS's
minimum spacing — which is what stops a fillet of many segments at a tiny
distance from emitting duplicate points.
=#
mutable struct _OffsetCurve{P}
    pts::Vector{P}
    min_spacing::Float64
end
_OffsetCurve(::Type{P}, d::Float64) where {P} = _OffsetCurve{P}(P[], 1e-4 * d)

function _offset_add_point!(c::_OffsetCurve{P}, p::P) where {P}
    if !isempty(c.pts)
        q = c.pts[end]
        hypot(p[1] - q[1], p[2] - q[2]) < c.min_spacing && return nothing
    end
    push!(c.pts, p)
    return nothing
end

function _offset_close_ring!(c::_OffsetCurve)
    isempty(c.pts) && return nothing
    c.pts[1] == c.pts[end] || push!(c.pts, c.pts[1])
    return nothing
end

#=
JTS `addDirectedFillet`: `n` points starting AT `a0` and sweeping toward `a1` in
direction `dir`, with the end point deliberately NOT emitted (the caller adds the
exact offset endpoint instead of a recomputed one). `n` is rounded, not ceiled,
so a sweep below half a quantum emits nothing at all.
=#
function _offset_directed_fillet!(c::_OffsetCurve{P}, centre::P, a0::Float64, a1::Float64,
        dir::Int, r::Float64, quantum::Float64) where {P}
    total = abs(a0 - a1)
    n = trunc(Int, total / quantum + 0.5)
    n < 1 && return nothing
    inc = total / n
    for i in 0:(n - 1)
        θ = a0 + dir * i * inc
        _offset_add_point!(c, (centre[1] + r * _offset_snap(cos(θ)),
                               centre[2] + r * _offset_snap(sin(θ))))
    end
    return nothing
end

#=
JTS `addCornerFillet`: the exact offset point `p0`, the arc, then the exact
offset point `p1`. The `±2π` normalization is what makes the sweep follow the
turn direction rather than the shorter way round.
=#
function _offset_corner_fillet!(c::_OffsetCurve{P}, centre::P, p0::P, p1::P,
        dir::Int, r::Float64, quantum::Float64) where {P}
    a0 = atan(p0[2] - centre[2], p0[1] - centre[1])
    a1 = atan(p1[2] - centre[2], p1[1] - centre[1])
    if dir == _OFFSET_CW
        a0 <= a1 && (a0 += 2π)
    else
        a0 >= a1 && (a0 -= 2π)
    end
    _offset_add_point!(c, p0)
    _offset_directed_fillet!(c, centre, a0, a1, dir, r, quantum)
    _offset_add_point!(c, p1)
    return nothing
end

# JTS `createCircle`: clockwise from `(x + d, y)`.
function _offset_circle(::Type{P}, p::P, d::Float64, quantum::Float64) where {P}
    c = _OffsetCurve(P, d)
    _offset_add_point!(c, (p[1] + d, p[2]))
    _offset_directed_fillet!(c, p, 0.0, 2π, _OFFSET_CW, d, quantum)
    _offset_close_ring!(c)
    return c.pts
end

# ## Input extraction and cleaning

function _offset_dedup_copy(pts::Vector{P}) where {P}
    out = P[]
    sizehint!(out, length(pts))
    for p in pts
        (isempty(out) || out[end] != p) && push!(out, p)
    end
    return out
end

_offset_points(::Type{P}, geom) where {P} =
    P[(Float64(GI.x(p)), Float64(GI.y(p))) for p in GI.getpoint(geom)]

# The distinct vertices of a ring, closing point removed.
function _offset_ring_vertices(::Type{P}, ring) where {P}
    v = _offset_dedup_copy(_offset_points(P, ring))
    length(v) > 1 && v[1] == v[end] && pop!(v)
    return v
end

#=
## Ring orientation

`_offset_ring_orientation` returns `+1` for a counter-clockwise ring, `-1` for a
clockwise one and `0` for one with no signed area. It is the only thing that
decides which side of a ring the offset goes on, so getting it wrong does not
perturb the result — it buffers the complement.

A plain shoelace sum over the raw coordinates does get it wrong. Every term is
`x_i·y_{i+1} - x_{i+1}·y_i`, a difference of two products of the *absolute*
coordinates, so the ring's area is computed as a difference of quantities of size
`coord²`. A unit square at `(1e8, 1e8)` has terms of size 1e16 and an area of 1:
the products round to the same Float64 and the sum comes out exactly zero, so the
shell is never reversed and `buffer(square, 0.1)` returned area 0.04 against
GEOS's 1.43. At `(1e10, 1e10)` both traversal orders fail.

So this is the same two-tier filter the polygon builder uses on arrangement rings
(`_ring_is_ccw_exact` in `maximal_edge_ring.jl`): translate to the first vertex,
which makes each term a product of *ring-sized* numbers, accumulate a forward
error bound alongside, and escalate to `Rational{BigInt}` — exact on
Float64-derived values — when the bound does not certify the sign. Translation is
exact whenever the ring spans less than the exponent range, and the bound covers
it when it is not.
=#
function _offset_ring_orientation(v::Vector{P}) where {P}
    length(v) < 3 && return 0
    (acc, bound) = _offset_ring_area2_bounded(v)
    abs(acc) > bound && return acc > 0 ? 1 : -1
    return sign(_offset_ring_area2_exact(v))
end

# The translated Float64 doubled signed area and a forward error bound on it:
# `n * mag` for the accumulated summation error, `cmax * len` for the ½-ulp
# rounding of each translated coordinate acting on the ring's L1 extent.
function _offset_ring_area2_bounded(v::Vector{P}) where {P}
    n = length(v)
    ox, oy = v[1]
    acc = 0.0; mag = 0.0; len = 0.0; cmax = max(abs(ox), abs(oy))
    ax = 0.0; ay = 0.0                       # translated previous vertex
    @inbounds for i in 2:n
        px, py = v[i]
        cmax = max(cmax, abs(px), abs(py))
        bx = px - ox; by = py - oy
        t1 = ax * by; t2 = ay * bx
        acc += t1 - t2
        mag += abs(t1) + abs(t2)
        len += abs(bx) + abs(by)
        ax = bx; ay = by
    end
    #-- the closing term (last -> first) has both translated factors zero
    u = 0.5 * eps(Float64)
    return (acc, 8 * u * (n * mag + cmax * len))
end

# The same sum in exact rational arithmetic. Float64 values convert to
# `Rational{BigInt}` exactly, so the sign returned here is the ring's true one.
function _offset_ring_area2_exact(v::Vector{P}) where {P}
    ox = Rational{BigInt}(v[1][1]); oy = Rational{BigInt}(v[1][2])
    acc = zero(Rational{BigInt})
    ax = zero(Rational{BigInt}); ay = zero(Rational{BigInt})
    @inbounds for i in 2:length(v)
        bx = Rational{BigInt}(v[i][1]) - ox
        by = Rational{BigInt}(v[i][2]) - oy
        acc += ax * by - ay * bx
        ax = bx; ay = by
    end
    return acc
end

#=
## The fully-eroded ring pre-check

JTS `BufferCurveSetBuilder.isRingFullyEroded`: a ring whose own erosion is
already empty contributes no curve at all. `_ring_fully_eroded` is that test on
open, deduplicated ring vertices, and `_buffer_parts` is where it is applied —
to a shell when `d < 0` and to a hole when `d > 0`, the only two sides that can
erode.

**Why the engine wants it even though it needs no collapse heuristics.** The
winding identity is unconditional, so the curve of a collapsed ring still gives
the right (empty) answer — it just costs Θ(n²) self-crossings to say so. A
1000-gon eroded past its apothem inverts into a tangle of 1.5 million crossings
and 69 s where the answer is `∅`.

**Why it is sound.** Every `p` in a ring's interior is within `envMin / 2` of the
ring: walking from `p` toward either end of the narrow envelope axis leaves the
ring's interior after at most that far, so it crosses the ring. Hence
`2|d| >= envMin` implies `dist(p, ring) <= |d|` for every interior point, i.e.
the ring's own erosion has empty interior. For a triangle the incircle radius is
the exact threshold, which is what `_triangle_fully_eroded` uses instead.

That is a statement about the *true* buffer, and it composes with the winding
identity `w = 1_P ± cov` (`buffer.jl`, PLAN.md §2) one ring at a time:

- **shell, `d < 0`** — the polygon's erosion is contained in its shell's erosion,
  which is empty, so the whole part vanishes: dropping the shell's curve *and*
  its holes' leaves `w ≡ 0` there, and `{w ≥ 1} = ∅` either way. Parts of a
  MultiPolygon or collection are independent: a dropped part's curves are
  supported on that part alone (its pieces lie inside it), so the winding field
  elsewhere is untouched.
- **hole, `d > 0`** — the hole's curve contributes `−1_hole + cov_hole`, and
  `cov_hole` is supported inside the closed hole, so dropping it changes `w` only
  there, from `cov_shell + cov_hole` to `cov_shell + 1`. The test says every
  point of the hole is within `|d|` of `∂hole ⊆ ∂P`, i.e. inside the dilation, so
  both fields are `≥ 1` and the kept region is the same. The shell's dilation
  fills the hole.

**It is sufficient, not necessary.** A ring can erode to nothing with a wide
envelope — a regular `n`-gon at exactly its apothem has `envMin = 2R > 2·apothem`
— and those cases go through the engine, which returns the right empty answer
after building the whole Θ(n²) tangle: 2.7 s for a 1000-gon.
=#
_ring_envelope_min_dim(v::Vector{P}) where {P} =
    min(maximum(p -> p[1], v) - minimum(p -> p[1], v),
        maximum(p -> p[2], v) - minimum(p -> p[2], v))

#=
JTS `isTriangleErodedCompletely`, a *precise* test: the inner buffer of a
triangle converges on its incentre, so it is empty exactly when the incircle
radius is below `|d|`. Written as JTS writes it — `Triangle.inCentre` weighted by
the opposite side lengths, then `Distance.pointToSegment` to one side — rather
than as `2·area / perimeter`, because the two disagree in the last ulp and this
is a test whose interesting inputs sit on its own boundary.
=#
function _triangle_fully_eroded(v::Vector{P}, d::Float64) where {P}
    a, b, c = v[1], v[2], v[3]
    l0 = hypot(b[1] - c[1], b[2] - c[2])
    l1 = hypot(a[1] - c[1], a[2] - c[2])
    l2 = hypot(a[1] - b[1], a[2] - b[2])
    circum = l0 + l1 + l2
    circum > 0 || return true
    cx = (l0 * a[1] + l1 * b[1] + l2 * c[1]) / circum
    cy = (l0 * a[2] + l1 * b[2] + l2 * c[2]) / circum
    return _offset_point_seg_dist(cx, cy, a, b) < abs(d)
end

# JTS `Distance.pointToSegment`.
function _offset_point_seg_dist(px::Float64, py::Float64, a::P, b::P) where {P}
    dx = b[1] - a[1]; dy = b[2] - a[2]
    len2 = dx * dx + dy * dy
    len2 == 0 && return hypot(px - a[1], py - a[2])
    r = ((px - a[1]) * dx + (py - a[2]) * dy) / len2
    r <= 0 && return hypot(px - a[1], py - a[2])
    r >= 1 && return hypot(px - b[1], py - b[2])
    return abs((a[2] - py) * dx - (a[1] - px) * dy) / sqrt(len2)
end

# `v` is open and deduplicated; the caller has already checked that this ring is
# on the eroding side of `d` (a shell with `d < 0`, a hole with `d > 0`).
function _ring_fully_eroded(v::Vector{P}, d::Float64) where {P}
    length(v) < 3 && return true
    length(v) == 3 && return _triangle_fully_eroded(v, d)
    #-- `>=`, where JTS has `>`: at `2|d| == envMin` every interior point is at
    #-- distance exactly `envMin / 2`, so the *open* erosion is still empty and the
    #-- closed one has no area. Measured to change no result across the 2247-case
    #-- collapse suite and the 28 harness workloads, and it is what takes an
    #-- n x n grid of unit squares eroded by exactly half a square off the
    #-- engine's slow path.
    return 2 * abs(d) >= _ring_envelope_min_dim(v)
end

# Shell CCW, holes CW: the material interior is then always on the LEFT of travel.
function _offset_canonical_polygon(::Type{P}, poly) where {P}
    shell = _offset_ring_vertices(P, GI.getexterior(poly))
    _offset_ring_orientation(shell) < 0 && reverse!(shell)
    holes = Vector{P}[]
    for h in GI.gethole(poly)
        hv = _offset_ring_vertices(P, h)
        _offset_ring_orientation(hv) > 0 && reverse!(hv)
        push!(holes, hv)
    end
    return shell, holes
end

# One flat (polygons, lines, points) decomposition of any geometry.
struct _BufferParts{P}
    polys::Vector{Tuple{Vector{P}, Vector{Vector{P}}}}
    lines::Vector{Vector{P}}
    points::Vector{P}
end
_BufferParts(::Type{P}) where {P} =
    _BufferParts{P}(Tuple{Vector{P}, Vector{Vector{P}}}[], Vector{P}[], P[])

#=
Abstract traits throughout, so `GI.Line`, `GeometryBasics.Triangle` and every
other alias lands in the right bucket without a method per spelling. The one
`AbstractGeometryCollectionTrait` method covers all four multi-geometries and
`GeometryCollection` alike, since the recursion only needs `GI.getgeom`.
=#
_offset_flatten!(b::_BufferParts{P}, g, ::GI.AbstractPointTrait) where {P} =
    push!(b.points, (Float64(GI.x(g)), Float64(GI.y(g))))
_offset_flatten!(b::_BufferParts{P}, g, ::GI.AbstractCurveTrait) where {P} =
    push!(b.lines, _offset_dedup_copy(_offset_points(P, g)))
_offset_flatten!(b::_BufferParts{P}, g, ::GI.AbstractPolygonTrait) where {P} =
    push!(b.polys, _offset_canonical_polygon(P, g))
function _offset_flatten!(b::_BufferParts, g, ::GI.AbstractGeometryCollectionTrait)
    for sub in GI.getgeom(g)
        _offset_flatten!(b, sub, GI.trait(sub))
    end
    return b
end
_offset_flatten!(b::_BufferParts, g, trait) = throw(ArgumentError(
    "buffer: no native buffer for a geometry with trait $(trait). Curved and " *
    "polyhedral-surface geometries are not supported; `buffer(GEOS(), geom, d)` " *
    "may be."))

#=
## Distances below the coordinate grid

Offsetting is `p ± d·u`, so a `|d|` below half an ulp of the largest coordinate
present cannot move any vertex off its own grid point: every offset point rounds
back onto the input, and what reaches the arrangement is a ring against a
bit-identical copy of itself — no interior, and no reason to build it. Such a
distance is therefore snapped to zero, where the answer is GEOS's: an areal
geometry comes back unchanged and anything of lower dimension is empty, since a
zero-radius sausage has no area.

Review found this boundary empirically — at Earth-scale coordinates a 1-ulp
buffer returned a valid polygon and a 0.1-ulp one returned the input. The snap
makes that a rule instead of a coincidence of where the rounding landed.

The scan is over the already-flattened parts, so it costs one pass over
coordinates that have just been read, and nothing at all once `d` is zero.
=#
function _offset_snap_distance(b::_BufferParts{P}, d::Float64) where {P}
    d == 0.0 && return d
    scale = 0.0
    for (shell, holes) in b.polys
        scale = _offset_max_coord(scale, shell)
        for h in holes
            scale = _offset_max_coord(scale, h)
        end
    end
    for l in b.lines
        scale = _offset_max_coord(scale, l)
    end
    scale = _offset_max_coord(scale, b.points)
    return abs(d) <= 0.5 * eps(scale) ? 0.0 : d
end

function _offset_max_coord(scale::Float64, v::Vector{P}) where {P}
    @inbounds for p in v
        scale = max(scale, abs(p[1]), abs(p[2]))
    end
    return scale
end

#=
Degenerate demotion, mirroring JTS `getRingCurve` / `getLineCurve`: a ring with
fewer than three distinct vertices is linework, and a line with one distinct
point is a point. For a negative distance only areal parts contribute — nothing
one-dimensional has an interior to erode — and degenerate shells drop out
entirely.
=#
function _buffer_parts(raw::_BufferParts{P}, d::Float64, prune::Bool) where {P}
    out = _BufferParts(P)
    if d < 0
        for (shell, holes) in raw.polys
            #-- a degenerate shell has no interior to erode, pre-check or not
            length(shell) < 3 && continue
            #-- a shell that erodes to nothing takes its holes with it
            prune && _ring_fully_eroded(shell, d) && continue
            push!(out.polys, (shell, [h for h in holes if length(h) >= 3]))
        end
        return out
    end
    lines = copy(raw.lines)
    for (shell, holes) in raw.polys
        if length(shell) < 3
            #-- a collapsed shell buffers like linework; its holes are then meaningless
            push!(lines, shell)
            continue
        end
        kept = Vector{P}[]
        for h in holes
            #-- a hole the dilation swallows needs no curve: the shell's fills it
            prune && length(h) >= 3 && d > 0 && _ring_fully_eroded(h, d) && continue
            length(h) < 3 ? push!(lines, h) : push!(kept, h)
        end
        push!(out.polys, (shell, kept))
    end
    append!(out.points, raw.points)
    for l in lines
        length(l) == 1 ? push!(out.points, l[1]) :
            length(l) >= 2 ? push!(out.lines, l) : nothing
    end
    return out
end

# ## The curve generator (JTS `OffsetCurveBuilder`, vertex-path variant)

#=
The state machine of JTS `OffsetSegmentGenerator`: three consecutive input
vertices `s0, s1, s2` and the offsets of the two segments they span.
=#
mutable struct _OffsetGen{P}
    curve::_OffsetCurve{P}
    d::Float64
    quantum::Float64
    side::Int
    s0::P
    s1::P
    s2::P
    off0::Tuple{P, P}
    off1::Tuple{P, P}
    #-- guard 3: `false` until this side has emitted one untrimmed vertex path
    armed::Bool
    #-- guard 2: squared distance from the current offset segment's start to the
    #-- last emitted point, when that point sits on the segment; 0 otherwise
    pos2::Float64
end
function _OffsetGen(::Type{P}, d::Float64, quantum::Float64) where {P}
    z = (0.0, 0.0)
    return _OffsetGen{P}(_OffsetCurve(P, d), d, quantum, _OFFSET_LEFT, z, z, z, (z, z), (z, z),
                         false, 0.0)
end

function _offset_init_side!(g::_OffsetGen{P}, s1::P, s2::P, side::Int) where {P}
    g.s1 = s1; g.s2 = s2; g.side = side
    g.off1 = _offset_segment(s1, s2, side, g.d)
    g.armed = false; g.pos2 = 0.0
    return nothing
end

@inline _offset_dist2(p::P, q::P) where {P} = (p[1] - q[1])^2 + (p[2] - q[2])^2

#=
The point where the offset segments `(a0, a1)` and `(b0, b1)` cross, or `nothing`
if they do not cross *properly*. "Properly" is JTS `addInsideTurn`'s segment test
tightened to a strict crossing: each segment's endpoints must be strictly on
opposite sides of the other's line, which is what puts the crossing in the
interior of both and so the removed kite inside both half-rectangles. The four
orientations are exact; the point itself is the plain parametric formula written
from `a1` (the offset endpoint at the corner) rather than from `a0`, so the
subtraction that produces it works on the short arm of the corner instead of on
two lengths of the input edge.
=#
function _offset_proper_crossing(a0::P, a1::P, b0::P, b1::P) where {P}
    o1 = rk_orient(Planar(), a0, a1, b0; exact = True())
    o2 = rk_orient(Planar(), a0, a1, b1; exact = True())
    ((o1 > 0) & (o2 < 0)) | ((o1 < 0) & (o2 > 0)) || return nothing
    o3 = rk_orient(Planar(), b0, b1, a0; exact = True())
    o4 = rk_orient(Planar(), b0, b1, a1; exact = True())
    ((o3 > 0) & (o4 < 0)) | ((o3 < 0) & (o4 > 0)) || return nothing
    ax = a1[1] - a0[1]; ay = a1[2] - a0[2]
    bx = b1[1] - b0[1]; by = b1[2] - b0[2]
    den = ax * by - ay * bx
    den == 0 && return nothing
    w = ((b0[1] - a1[1]) * by - (b0[2] - a1[2]) * bx) / den
    x = a1[1] + w * ax; y = a1[2] + w * ay
    (isfinite(x) & isfinite(y)) || return nothing
    return (x, y)
end

_offset_add_last!(g::_OffsetGen) = _offset_add_point!(g.curve, g.off1[2])

function _offset_add_next!(g::_OffsetGen{P}, p::P) where {P}
    g.s0 = g.s1; g.s1 = g.s2; g.s2 = p
    g.off0 = _offset_segment(g.s0, g.s1, g.side, g.d)
    g.off1 = _offset_segment(g.s1, g.s2, g.side, g.d)
    prev2 = g.pos2
    #-- every branch but a trim leaves the curve at the new segment's own start
    g.pos2 = 0.0
    g.s1 == g.s2 && return nothing
    o = Int(sign(rk_orient(Planar(), g.s0, g.s1, g.s2; exact = True())))
    if o == 0
        dot = (g.s1[1] - g.s0[1]) * (g.s2[1] - g.s1[1]) +
              (g.s1[2] - g.s0[2]) * (g.s2[2] - g.s1[2])
        #-- straight continuation: the two offsets are collinear, so nothing is needed
        dot > 0 && return nothing
        #-- exact reversal: a 180° fillet around the tip, away from the offset side
        _offset_corner_fillet!(g.curve, g.s1, g.off0[2], g.off1[1],
                               g.side == _OFFSET_LEFT ? _OFFSET_CW : _OFFSET_CCW,
                               g.d, g.quantum)
    elseif (o < 0 && g.side == _OFFSET_LEFT) || (o > 0 && g.side == _OFFSET_RIGHT)
        if _offset_dist2(g.off0[2], g.off1[1]) < (_OFFSET_SEPARATION_FACTOR * g.d)^2
            #-- near-parallel outside turn: one point, from the longer segment,
            #-- which is the one whose offset the join displaces least
            _offset_add_point!(g.curve, _offset_dist2(g.s0, g.s1) > _offset_dist2(g.s1, g.s2) ?
                                        g.off0[2] : g.off1[1])
        else
            #-- outside turn: fillet sweeping in the turn direction
            _offset_corner_fillet!(g.curve, g.s1, g.off0[2], g.off1[1], o, g.d, g.quantum)
        end
    else
        #-- inside turn: the single crossing point of the two offset segments when
        #-- all three guards in the header hold, else the path through the vertex.
        x = g.armed ? _offset_proper_crossing(g.off0[1], g.off0[2], g.off1[1], g.off1[2]) :
            nothing
        if x !== nothing && _offset_dist2(x, g.off0[1]) > prev2
            _offset_add_point!(g.curve, x)
            g.pos2 = _offset_dist2(x, g.off1[1])
        else
            g.armed = true
            _offset_add_point!(g.curve, g.off0[2])
            _offset_add_point!(g.curve, g.s1)
            _offset_add_point!(g.curve, g.off1[1])
        end
    end
    return nothing
end

# JTS `addLineEndCap`, round only: left offset end, clockwise half-circle, right end.
function _offset_end_cap!(g::_OffsetGen{P}, p0::P, p1::P) where {P}
    oleft = _offset_segment(p0, p1, _OFFSET_LEFT, g.d)
    oright = _offset_segment(p0, p1, _OFFSET_RIGHT, g.d)
    ang = atan(p1[2] - p0[2], p1[1] - p0[1])
    _offset_add_point!(g.curve, oleft[2])
    _offset_directed_fillet!(g.curve, p1, ang + π / 2, ang - π / 2, _OFFSET_CW, g.d, g.quantum)
    _offset_add_point!(g.curve, oright[2])
    return nothing
end

#=
JTS `computeLineBufferCurve`: the left side forward, an end cap, the left side of
the reversed line, a start cap. That traversal is clockwise around the buffer
region, so the ring is reversed once at the end to make it positive.
=#
function _offset_line_curve(::Type{P}, v::Vector{P}, d::Float64, quantum::Float64) where {P}
    g = _OffsetGen(P, d, quantum)
    n = length(v)
    sizehint!(g.curve.pts, 4n + 32)   # two sides, ~2 points per vertex, plus caps
    _offset_init_side!(g, v[1], v[2], _OFFSET_LEFT)
    for i in 3:n
        _offset_add_next!(g, v[i])
    end
    _offset_add_last!(g)
    _offset_end_cap!(g, v[n - 1], v[n])
    _offset_init_side!(g, v[n], v[n - 1], _OFFSET_LEFT)
    for i in (n - 2):-1:1
        _offset_add_next!(g, v[i])
    end
    _offset_add_last!(g)
    _offset_end_cap!(g, v[2], v[1])
    _offset_close_ring!(g.curve)
    return reverse(g.curve.pts)
end

#=
JTS `computeRingBufferCurve`, travelling in the ring's canonical direction.

One deviation, and it is structural: JTS passes `addStartPoint = false` for the
first corner so the closing segment can skip `off0.p1`. With the vertex path that
would replace a radial edge by a diagonal and break the winding identity, so the
start point is always emitted.
=#
function _offset_ring_curve(::Type{P}, v::Vector{P}, side::Int, d::Float64,
        quantum::Float64) where {P}
    g = _OffsetGen(P, d, quantum)
    n = length(v)
    sizehint!(g.curve.pts, 2n + 16)   # ~2 points per vertex when joins are short
    _offset_init_side!(g, v[n], v[1], side)
    for i in 2:n
        _offset_add_next!(g, v[i])
    end
    _offset_add_next!(g, v[1])
    _offset_close_ring!(g.curve)
    return g.curve.pts
end

#=
Every raw offset curve of `geom` at signed distance `d`, as closed rings oriented
positively around the buffer region.

Polygons are canonicalised so their interior is on the LEFT of travel, which
fixes the offset side from the sign of `d` alone: RIGHT for `d > 0` (away from the
interior), LEFT for `d < 0` (into it). The curve then keeps the ring's own travel
direction and needs no reversal. The polygon's own rings are NOT emitted — the
curve alone carries the whole region, verified across the 536-case corpus
including total erosion and a hole swallowing its shell.

`prune_eroded_rings` gates the fully-eroded pre-check only; the curves are the
same either way for every ring the check does not drop.
=#
function _raw_offset_curves(::Type{P}, geom, d::Float64, quadsegs::Int;
        prune_eroded_rings::Bool = true) where {P}
    flat = _BufferParts(P)
    _offset_flatten!(flat, geom, GI.trait(geom))
    d = _offset_snap_distance(flat, d)
    quantum = (π / 2) / max(quadsegs, 1)
    ad = abs(d)
    side = d > 0 ? _OFFSET_RIGHT : _OFFSET_LEFT
    parts = _buffer_parts(flat, d, prune_eroded_rings)
    rings = Vector{P}[]
    for (shell, holes) in parts.polys
        push!(rings, _offset_ring_curve(P, shell, side, ad, quantum))
        for h in holes
            push!(rings, _offset_ring_curve(P, h, side, ad, quantum))
        end
    end
    for l in parts.lines
        push!(rings, _offset_line_curve(P, l, ad, quantum))
    end
    for p in parts.points
        push!(rings, reverse(_offset_circle(P, p, ad, quantum)))
    end
    return rings
end
