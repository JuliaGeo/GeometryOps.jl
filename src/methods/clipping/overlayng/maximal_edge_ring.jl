# # Result ring linking — maximal rings, minimal-ring split, edge rings
#
# Phase 2b of the OverlayNG port (design doc §3). Ports two JTS files kept
# together because they form one pipeline over the marked result-area edges:
#   - `MaximalEdgeRing.java`  → maximal-ring linking + the minimal-ring split
#     at self-touching nodes (`_MaxEdgeRing` + `_link_result_area_max_ring_at_node!`).
#   - `OverlayEdgeRing.java`  → one minimal result ring: its coordinates, its
#     shell/hole role (decided by exact kernel predicates over the arrangement's
#     nodes — `_ring_is_ccw_exact`), and hole containment (design §3
#     amendment 5 — never planar even-odd on emitted coordinates; see the KNOWN
#     GAP note there, containment is still decided on emitted coordinates).
#
# Ring linkage lives in the phase-2a `OverlayEdge` handles: `next_result_max`
# links maximal rings, `next_result` links minimal rings, and `max_edge_ring` /
# `edge_ring` are integer handles (`0` = null) into a builder's ring vectors.
# Each `NodedEdge`/`OverlayEdge` is a single straight segment between two nodes,
# so a ring's coordinates are exactly the sequence of its edges' node points —
# no intermediate `addCoordinates` bookkeeping is needed.
#
# Everything here is internal to GeometryOps — nothing is exported.

# A maximal edge ring: a cycle of result-area half-edges linked by
# `next_result_max`. Identified by `id` (stored in each member's `max_edge_ring`
# handle) so `_attach_max_edges!` can detect revisits (port of JTS
# `MaximalEdgeRing`).
mutable struct _MaxEdgeRing
    id::Int32
    start_edge::Int32
end

# One minimal result ring (port of JTS `OverlayEdgeRing`): the emitted output
# coordinates (`ring_pts`, closed), the arrangement node ids the ring visits
# (`node_ids`, open — the ring's EXACT identity, from which the shell/hole role
# is decided), its shell/hole role, a bounding box for extent pruning, its
# assigned shell / contained holes (handles, `0` = null), and a lazily-built
# indexed point-in-area locator over its own ring.
#
# `node_ids` and `ring_pts` are not in bijection: two distinct nodes can realize
# to the same output coordinate, which `_ring_add!` drops from `ring_pts` (see
# below). That is exactly why the role is taken from `node_ids`.
mutable struct _OverlayEdgeRing{P}
    id::Int32
    start_edge::Int32
    ring_pts::Vector{Tuple{Float64, Float64}}
    node_ids::Vector{Int32}
    is_hole::Bool
    bbox::NTuple{4, Float64}   # (xmin, xmax, ymin, ymax) of ring_pts
    shell::Int32
    holes::Vector{Int32}
    locator::Any               # Union{Nothing, IndexedPointInAreaLocator}, lazy
end

# The polygon-builder working context (JTS `PolygonBuilder`'s mutable state). Held
# together so the ring-linking and placement functions share the graph edge store,
# the arrangement (for `node_point`), the manifold/exact predicate context, and the
# growing ring collections. Parameterized on `M`/`P`/`E` so `m`/`exact` stay
# concrete and the `Planar`/`Spherical` methods dispatch.
mutable struct _PolyBuilderCtx{M <: Manifold, P, E}
    m::M
    edges::Vector{OverlayEdge{P}}
    arr::NodedArrangement{P}
    exact::E
    max_rings::Vector{_MaxEdgeRing}
    edge_rings::Vector{_OverlayEdgeRing{P}}
    shell_list::Vector{Int32}      # handles into edge_rings
    free_hole_list::Vector{Int32}
    #-- free holes whose OWN shell was dropped as sub-grid, so they have no
    #-- containing shell to be entitled to. Kept as a separate list rather than a
    #-- flag on the ring because it records why the hole is free, which is what
    #-- decides whether being unplaceable is a defect or the expected outcome
    #-- (`_place_free_holes!`).
    orphan_hole_list::Vector{Int32}
    #-- per-node memo of the exact kernel position used by the orientation
    #-- predicate (spherical only; see `_node_kernel_point`). Allocated on
    #-- first use so the planar path pays nothing.
    kernel_cache::Vector{P}
    kernel_ok::Vector{Bool}
end

_PolyBuilderCtx(m::M, edges::Vector{OverlayEdge{P}}, arr::NodedArrangement{P}, exact::E,
        max_rings, edge_rings, shell_list, free_hole_list) where {M <: Manifold, P, E} =
    _PolyBuilderCtx{M, P, E}(m, edges, arr, exact, max_rings, edge_rings, shell_list,
                             free_hole_list, Int32[], P[], Bool[])

@inline _ctx_point_type(::_PolyBuilderCtx{M, P}) where {M, P} = P

# ## Maximal-ring linking at a node (port of `linkResultAreaMaxRingAtNode`)
#
# Design §3 amendment 4 (the known trap): this is called UNGATED for every
# in-result edge (JTS's `// TODO: skip already-linked` is deliberately
# unfulfilled — gating on already-linked loses degree-2 nodes whose lone
# out-edge was pre-linked as an in-edge). The per-node scan's own early return
# on an already-linked in-edge provides the necessary idempotency.
function _link_result_area_max_ring_at_node!(edges, node_edge::Integer)
    #-- precondition: node_edge is in the result area
    end_out = he_onext(edges, node_edge)
    curr_out = end_out
    #-- state machine: 1 = find an incoming result edge, 2 = link to an outgoing
    state = 1
    curr_result_in = Int32(0)
    while true
        #-- if the found in-edge is already linked, this node is done
        (curr_result_in != 0 && oe_is_result_max_linked(edges, curr_result_in)) && return nothing
        if state == 1
            curr_in = he_sym(edges, curr_out)
            if oe_in_result_area(edges, curr_in)
                curr_result_in = curr_in
                state = 2
            end
        else # state == 2
            if oe_in_result_area(edges, curr_out)
                oe_set_next_result_max!(edges, curr_result_in, curr_out)
                state = 1
            end
        end
        curr_out = he_onext(edges, curr_out)
        curr_out == end_out && break
    end
    state == 2 && throw(_OverlayTopologyError("no outgoing edge found"))
    return nothing
end

# Attach the edges of a maximal ring, tagging each with `mr.id` (port of the
# `MaximalEdgeRing` constructor / `attachEdges`).
function _attach_max_edges!(ctx, mr::_MaxEdgeRing)
    edges = ctx.edges
    edge = mr.start_edge
    while true
        edge == 0 && throw(_OverlayTopologyError("Ring edge is null"))
        edges[edge].max_edge_ring == mr.id &&
            throw(_OverlayTopologyError("Ring edge visited twice at max-ring build"))
        oe_next_result_max(edges, edge) == 0 &&
            throw(_OverlayTopologyError("Ring edge missing at max-ring build"))
        edges[edge].max_edge_ring = mr.id
        edge = oe_next_result_max(edges, edge)
        edge == mr.start_edge && break
    end
    return nothing
end

# ## Minimal-ring split (ports of `buildMinimalRings` / `linkMinimalRings`)
#
# Splits a self-touching maximal ring into OGC-valid minimal rings by relinking
# the max-ring edges in the OPPOSITE (CW) orientation via `next_result`. This is
# exactly the piece the spike prototypes faked (which produced invalid unions on
# many-island geometries).

# Build the minimal rings of a maximal ring, returning their handles.
function _build_minimal_rings!(ctx, mr::_MaxEdgeRing)
    _link_minimal_rings!(ctx, mr)
    min_rings = Int32[]
    edges = ctx.edges
    e = mr.start_edge
    while true
        edges[e].edge_ring == 0 && push!(min_rings, _new_edge_ring!(ctx, e))
        e = oe_next_result_max(edges, e)
        e == mr.start_edge && break
    end
    return min_rings
end

function _link_minimal_rings!(ctx, mr::_MaxEdgeRing)
    edges = ctx.edges
    e = mr.start_edge
    while true
        _link_min_ring_edges_at_node!(ctx, e, mr)
        e = oe_next_result_max(edges, e)
        e == mr.start_edge && break
    end
    return nothing
end

# Port of `linkMinRingEdgesAtNode`: relink this max ring's edges around one node
# into minimal rings (CW orientation, via `next_result`).
function _link_min_ring_edges_at_node!(ctx, node_edge::Integer, mr::_MaxEdgeRing)
    edges = ctx.edges
    end_out = Int32(node_edge)
    curr_max_ring_out = Int32(node_edge)
    curr_out = he_onext(edges, node_edge)
    while true
        _is_already_linked_min(edges, he_sym(edges, curr_out), mr) && return nothing
        if curr_max_ring_out == 0
            curr_max_ring_out = _select_max_out_edge(edges, curr_out, mr)
        else
            curr_max_ring_out = _link_max_in_edge!(edges, curr_out, curr_max_ring_out, mr)
        end
        curr_out = he_onext(edges, curr_out)
        curr_out == end_out && break
    end
    curr_max_ring_out != 0 &&
        throw(_OverlayTopologyError("Unmatched edge found during min-ring linking"))
    return nothing
end

@inline _is_already_linked_min(edges, edge::Integer, mr::_MaxEdgeRing) =
    edges[edge].max_edge_ring == mr.id && oe_is_result_linked(edges, edge)

@inline _select_max_out_edge(edges, curr_out::Integer, mr::_MaxEdgeRing) =
    edges[curr_out].max_edge_ring == mr.id ? Int32(curr_out) : Int32(0)

@inline function _link_max_in_edge!(edges, curr_out::Integer, curr_max_ring_out::Integer,
        mr::_MaxEdgeRing)
    curr_in = he_sym(edges, curr_out)
    edges[curr_in].max_edge_ring != mr.id && return Int32(curr_max_ring_out)
    oe_set_next_result!(edges, curr_in, curr_max_ring_out)
    return Int32(0)
end

# ## Minimal ring construction (port of the `OverlayEdgeRing` constructor)

#=
Port of JTS `CoordinateList.add(coord, allowRepeated = false)`, the accumulator
`OverlayEdge.addCoordinates` fills a ring's coordinates through: a point equal
to the current last one is not appended.

This is not cosmetic here. The arrangement's nodes are symbolic and exact, so two
*distinct* nodes can realize to the same Float64 coordinate at emission (design
§2.6 rounds once, there); the ring then carries a repeated point, which is not a
legal `LinearRing` vertex sequence. JTS never sees the situation because its
noder rounds while noding, so the two nodes are already one.
=#
@inline function _ring_add!(pts::Vector{Tuple{Float64, Float64}}, p::Tuple{Float64, Float64})
    (isempty(pts) || pts[end] != p) && push!(pts, p)
    return nothing
end

#=
## Ring collapse at emission

A minimal ring is *collapsed* when no Float64 image of it is a legal, faithful
`LinearRing`. That is the emission-side analogue of JTS's `DIM_COLLAPSE` edge
label, arrived at one stage later because this engine rounds at emission rather
than while noding (design §2.6). Collapsed rings are excluded from the result for
the reason JTS excludes a collapse: they contribute no representable area, and
keeping them produces an invalid result geometry.

There are two halves, and neither implies the other:

1. `_ring_image_is_degenerate` — the emitted image has fewer than three distinct
   vertices, so it is not a legal `LinearRing` at all (JTS's
   `GeometryFactory.createLinearRing` rejects it outright). This is the ring
   whose two sides are exactly distinct but round to the same output
   coordinates. Exact equality of emitted coordinates; no tolerance.

2. `_ring_is_subgrid` — the ring is exactly non-degenerate AND its emitted image
   has three or more distinct vertices, yet its exact mean width is finer than
   the output format's grid step where the ring sits. Its image is then a ring
   whose sides have been displaced past each other by rounding: a self-touching
   or self-crossing needle, invalid output from a correct exact face.

Of the 7 spike holes that survived the spherical watershed cascade before half 2
existed, half 1 catches ZERO — every one of them emits three or more distinct
vertices, and two emit a perfectly ordinary-looking ring with a positive area
that `LG.isValid` accepts. Half 2 catches 5 of the 7 outright; the remaining two
are 30 µm-wide triangles that never form once the other five stop being emitted
into the next level of the cascade, so the end-to-end result is 0 of 7.

That does not make half 1 redundant — it is the only one of the two that needs no
arithmetic, and a ring with two distinct vertices has no legal image whatever its
width — but it does mean half 1 was never going to be enough. It fires on none of
the 2974 minimal rings the watershed cascade builds.

`_ring_is_collapsed(ctx, r)` is the decision the builder makes. The
one-argument method is half 1 only — a property of the ring alone, kept because
it is exactly what a caller holding no builder context can answer.
=#
@inline _ring_image_is_degenerate(r::_OverlayEdgeRing) = length(r.ring_pts) < 4

@inline _ring_is_collapsed(r::_OverlayEdgeRing) = _ring_image_is_degenerate(r)

@inline _ring_is_collapsed(ctx, r::_OverlayEdgeRing) =
    _ring_image_is_degenerate(r) || _ring_is_subgrid(ctx, r)

#=
### The sub-grid half, and the one magnitude-relative threshold in the engine

Design §0 says no decision reads a constructed coordinate, and the engine carries
no tolerance. This test is the one place a length scale appears, and it is worth
being precise about why it is not a tolerance in the sense §0 forbids.

A tolerance is a claim about the GEOMETRY: "features closer than ε are the same
feature". That claim is what snapping engines make, and it is what this port
refuses — it is why the arrangement is exact and why `_ring_is_ccw_exact` decides
on node ids rather than on `ring_pts`.

`u` below is instead a property of the OUTPUT FORMAT. `node_point` emits
`Tuple{Float64,Float64}`, and at coordinate magnitude `m` that format's
representable positions are `eps(m)` apart. A ring whose exact width is a
fraction of `eps(m)` has no faithful image in the format: every candidate image
either merges its two sides (half 1) or pushes them past one another (an invalid
ring). Dropping it is not deciding that the ring is "too small to be real" — the
exact arrangement is not consulted about whether it is real, and it IS real. It
is deciding that the result geometry cannot say so, which is a statement about
`Vector{Tuple{Float64,Float64}}`, not about the sphere or the plane.

Consequently `u` is derived from where the ring sits, never from a user tolerance
or a global constant, and it shrinks to nothing near the coordinate origin —
exactly where the format really can resolve arbitrarily fine detail, and exactly
where this test correspondingly declines to fire. That is the signature of a
format property rather than a tolerance.

The decision itself still runs on node ids and kernel points (design §0): `u`
fixes only the yardstick, and `ring_pts` is read only to find the ring's
coordinate magnitude — never to measure the ring.

THE TEST. Drop the ring when its exact mean width is below `_RING_GRID_MARGIN`
grid steps. Writing `D` for the larger side of the ring's bounding box and `A`
for its exact area, that is

    |A_exact| < _RING_GRID_MARGIN · u · D

`A/D` is the ring's area spread over its own length — its mean thickness across
the long axis. (It stands in for the mean width `2A/P`, which the two are equal
to for a needle since such a curve runs out and back, `P = 2D`. `P` itself is
deliberately not used; see `_ring_width_below` for why it cannot be.) The margin
is 4 because emission displaces each vertex by up to ½ ulp per coordinate, so two
facing vertices can move up to `√2 · u ≈ 1.41 u` relative to each other; 4 is the
smallest power of two strictly above that.

CALIBRATION, and an honest account of what it does and does not show. Measured on
the Vancouver watersheds (1384 polygons, cascaded union over an STR leaf order),
sweeping the margin:

    margin   0.0    0.5    1.0    1.5    2.0    2.5    3.0    4.0
    holes      7      5      2      1      1      1      0      0
    valid  FALSE  FALSE   true   true   true   true   true   true

Two different bars, and they are not in the same place:

  * the CONTRACT — never emit invalid output — is met from margin 1 up;
  * agreement with the planar answer's census (1 polygon, 0 holes) needs 3.

4 is the derived value above, and it clears the empirical requirement of 3 with
33% to spare. That headroom is the whole reason to prefer the derivation to the
measurement: the measurement is one corpus.

WHAT IT COSTS, because it is not free. Faces narrower than the output grid are
not rare noise that only broken inputs produce — GEOS emits them routinely, and
dropping them means we no longer match GEOS coordinate for coordinate:

  * planar real data (watersheds, 4 corpus/order combinations, ~12000 minimal
    rings): ZERO rings dropped. Real projected data does not come near this.
  * the 1600-op fuzz sweep: 42 ops (2.6%) now differ from GEOS. Every one stays
    valid, every one is inside the suite's own rounding band, none diverges —
    but they are differences, and their widths run from 0.44 to 3.98 grid steps,
    i.e. straight through the range this test drops.
  * JTS `TestOverlayMisc` case 2: we used to equal GEOS exactly and now differ by
    one splinter of area 1.8e-22 — 0.0014 long and 2.7e-19 wide, six thousand
    times narrower than the grid it is written on.

So there is NO empty band between "spike" and "legitimate": the two populations
overlap, and every margin above zero trades GEOS parity on sub-grid splinters for
validity on sub-grid needles. This constant is a policy choice about which of the
two matters more, argued from the output format — not a separation discovered in
the data. Anything that reports it as the latter is overclaiming.
=#
const _RING_GRID_MARGIN = 4.0

#=
The output-grid step at the ring's coordinate magnitude, in the units the exact
width is measured in.

On both manifolds the grid is ANISOTROPIC — the step differs between the two
axes — and both take the SMALLER of the two steps. That is the conservative
reading: "is the ring narrower than the finest detail the format can resolve
anywhere here". It under-fires (a needle aligned with the coarser axis is
missed) rather than over-fires, and over-firing is the direction that would lose
real geometry.

Planar: `node_point` emits the coordinate itself, so the step is `eps` at the
ring's largest magnitude on each axis, and widths are already in those units.

Spherical: `node_point` emits (lon, lat) in DEGREES while the exact width is
measured on the unit sphere in radians. One `eps(lat)` step of latitude is
`eps(lat)·π/180` radians of arc; one `eps(lon)` step of longitude is
`eps(lon)·cos(φ)·π/180`.

Reading `ring_pts` here is deliberate and is not a §0 violation: the question
asked of it is "how far apart are representable Float64s in this neighbourhood",
which only the emitted coordinates can answer.
=#
_ring_grid_step(ctx::_PolyBuilderCtx, pts::Vector{Tuple{Float64, Float64}}) =
    _ring_grid_step(ctx.m, pts)

function _ring_grid_step(::Planar, pts::Vector{Tuple{Float64, Float64}})
    xm = 0.0; ym = 0.0
    @inbounds for p in pts
        xm = max(xm, abs(p[1])); ym = max(ym, abs(p[2]))
    end
    return min(eps(xm), eps(ym))
end

function _ring_grid_step(::Spherical, pts::Vector{Tuple{Float64, Float64}})
    lonmax = 0.0; latmax = 0.0; latsum = 0.0
    @inbounds for p in pts
        lonmax = max(lonmax, abs(p[1])); latmax = max(latmax, abs(p[2]))
        latsum += p[2]
    end
    latmid = latsum / length(pts)
    return min(eps(lonmax) * abs(cosd(latmid)), eps(latmax)) * (π / 180)
end

function _ring_is_subgrid(ctx, r::_OverlayEdgeRing)
    ids = r.node_ids
    length(ids) < 3 && return true
    u = _ring_grid_step(ctx, r.ring_pts)
    (u > 0 && isfinite(u)) || return false
    return _ring_width_below(ctx, ids, _RING_GRID_MARGIN * u)
end

#=
Whether `2·|A_exact| < thr · P` over the ring's exact node positions, in the
filter/escalate shape every other predicate in the port uses
(`_ring_is_ccw_exact`, `rk_orient`, …): a certified Float64 pass first, exact
arithmetic only for what it cannot decide.

Plain Float64 CANNOT decide these rings on its own and must not be trusted to.
The areas in question are ~1e-22 while the terms summing to them are ~1e-5 — 17
digits of cancellation. (Measured the hard way: a first calibration pass that
computed spherical ring areas in Float64 reported ratios two to four orders of
magnitude off, pure rounding noise, and showed no separation at all between the
spike population and the legitimate one. The separation above is only visible in
extended precision.) So the escalation is not an unlikely fallback here; it is
where every ring near the threshold is actually decided. It is affordable because
the filter sends only near-threshold rings there, and those have 3–4 nodes.

P DOES NOT APPEAR. The perimeter is a sum of square roots: it has no rational
form to escalate to, and its Float64 error is driven by the ½-ulp displacement of
each endpoint, i.e. by the ring's COORDINATE MAGNITUDE rather than by its size.
That is fatal for exactly the rings this test exists for. A three-node ring
spanning 3 ulps at coordinates of magnitude 3.3e6 has a perimeter of ~4e-9 and an
endpoint-driven error bound of ~9e-9 — the bound exceeds the quantity, no positive
lower bound on `P` is certifiable at all, and a version of this test that needed
one simply gave up and kept the ring. (That is not hypothetical: it kept a
1.1e-19-area ring in JTS `TestOverlay-misc-4` case 5, which then orphaned itself
and turned an invalid result into a THROW.)

So `P` is replaced by the tightest lower bound with a rational form. `D` is the
larger side of the ring's bounding box; any closed curve traverses each of its
extents twice, so

    P ≥ 2·D           and the test becomes    |A2_exact| < 2 · thr · D

writing `A2 = 2A` for the doubled area the port already computes. For a needle —
which is what fires — the curve runs out and back along its length, so `P = 2D` to
first order and the substitution is not a substitution at all; for a compact ring
it is up to ~1.4x conservative, in the safe direction. `D` is O(n) to compute,
exact as a rational, and needs no square root: the escalated comparison is

    A2²  <  4 · thr² · D²

Both branches remain one-sided and both err towards KEEPING the ring: the filter
decides only when `|acc| ± bound` clears the threshold outright, and the escalation
is exact. A ring is dropped only when it is *certified* sub-grid.
=#
function _ring_width_below(ctx::_PolyBuilderCtx{<:Planar}, ids::Vector{Int32}, thr::Float64)
    arr = ctx.arr
    (per, per_err) = _ring_perimeter_filter(arr, ids)
    (acc, bound) = _ring_area2_bounded(arr, ids)
    #-- `2|A| = |area2|`, and the true value lies in `|acc| ± bound`
    if per > per_err
        abs(acc) - bound > thr * (per + per_err) && return false
        abs(acc) + bound < thr * (per - per_err) && return true
    end
    #-- escalate: exact area against a certified rational lower bound on P
    a2 = _ring_area2_exact(arr, ids)
    plo = _ring_perimeter_lo_exact(arr, ids)
    return abs(a2) < Rational{BigInt}(thr) * plo
end

#=
The ring's perimeter, twice: a cheap Float64 estimate with an error bound, and a
certified lower bound on the exact one. Two routines rather than one because they
fail on opposite inputs.

The FILTER is one `sqrt` per edge over `node_point`, but its error is driven by
the ½-ulp displacement of each endpoint — by the ring's COORDINATE MAGNITUDE
rather than by its size. On a three-node ring spanning 3 ulps at coordinates of
magnitude 3.3e6 the perimeter is ~4e-9 and the bound is ~9e-9: the bound exceeds
the quantity and nothing can be certified either way. That is not a corner case,
it is precisely the shape this test exists to catch. (Measured: a version that
gave up there kept a 1.1e-19-area ring in JTS `TestOverlay-misc-4` case 5, which
then orphaned itself and turned an invalid result into a THROW.)

So when the filter cannot certify, the ESCALATION rebuilds the perimeter from the
arrangement's EXACT node coordinates, where there is no endpoint displacement to
account for at all. The exact perimeter is irrational — a sum of square roots —
but a certified rational LOWER bound suffices, because `P` sits on the larger side
of the comparison and under-stating it can only make the test fire less often.
Each root is taken in Float64 and backed off two ulps, covering both the
rational-to-Float64 conversion and the `sqrt`.

(A version that replaced `P` with the root-free closed form
`2·max(bbox extent, longest edge)` was tried and abandoned. Both terms are
genuine lower bounds on `P`, and the form is exactly tight for an unsubdivided
axis-aligned needle — but loose for anything else, measured at 0.76 on the
spherical synthetic, whose needle is both diagonal in R³ and subdivided at its
midpoint. It moved the drop boundary apart between the two manifolds. A threshold
that is meant to be a property of the output grid must not depend on which way the
ring points or how many vertices sit along its side.)
=#
function _ring_perimeter_filter(arr, ids::Vector{Int32})
    n = length(ids)
    (x0, y0) = node_point(arr, ids[1])
    px = x0; py = y0
    per = 0.0; cmax = max(abs(x0), abs(y0))
    @inbounds for i in 2:n
        (qx, qy) = node_point(arr, ids[i])
        cmax = max(cmax, abs(qx), abs(qy))
        dx = qx - px; dy = qy - py
        per += sqrt(dx * dx + dy * dy)
        px = qx; py = qy
    end
    dx = x0 - px; dy = y0 - py
    per += sqrt(dx * dx + dy * dy)
    return (per, 4 * eps(Float64) * n * (cmax + per))
end

# A certified rational lower bound on `sqrt(d2)`.
function _sqrt_lo(d2::Rational{BigInt})
    d2 <= 0 && return zero(Rational{BigInt})
    s = sqrt(Float64(d2))
    (isfinite(s) && s > 0) || return zero(Rational{BigInt})
    lo = prevfloat(s, 2)
    lo > 0 || return zero(Rational{BigInt})
    return Rational{BigInt}(lo)
end

function _ring_perimeter_lo_exact(arr, ids::Vector{Int32})
    n = length(ids)
    (x0, y0) = _exact_node_point(arr.nodes.keys[ids[1]])
    px = x0; py = y0
    acc = zero(Rational{BigInt})
    @inbounds for i in 2:n
        (qx, qy) = _exact_node_point(arr.nodes.keys[ids[i]])
        dx = qx - px; dy = qy - py
        acc += _sqrt_lo(dx * dx + dy * dy)
        px = qx; py = qy
    end
    dx = x0 - px; dy = y0 - py
    return acc + _sqrt_lo(dx * dx + dy * dy)
end

# Float64 doubled signed area of the emitted ring, with the certificate above.
_ring_area2_filter(arr, ids::Vector{Int32}) =
    ((acc, bound) = _ring_area2_bounded(arr, ids); (acc, abs(acc) > bound))

function _new_edge_ring!(ctx, start::Integer)
    P = _ctx_point_type(ctx)
    id = Int32(length(ctx.edge_rings) + 1)
    ring = _OverlayEdgeRing{P}(id, Int32(start), Tuple{Float64, Float64}[], Int32[],
                               false, (0.0, 0.0, 0.0, 0.0), Int32(0), Int32[], nothing)
    push!(ctx.edge_rings, ring)
    _compute_ring!(ctx, ring)
    return id
end

# Port of `computeRingPts` + `computeRing`: walk the minimal ring via
# `next_result`, collecting the arrangement node ids it visits and their emitted
# points; then derive the shell/hole role and the bounding box.
function _compute_ring!(ctx, ring::_OverlayEdgeRing)
    edges = ctx.edges
    pts = Tuple{Float64, Float64}[]
    ids = Int32[]
    origin = he_origin(edges, ring.start_edge)
    push!(ids, Int32(origin))
    _ring_add!(pts, node_point(ctx.arr, origin))
    edge = ring.start_edge
    while true
        edges[edge].edge_ring == ring.id &&
            throw(_OverlayTopologyError("Edge visited twice during ring-building"))
        dest = he_dest(edges, edge)
        edges[edge].edge_ring = ring.id
        ne = oe_next_result(edges, edge)
        ne == 0 && throw(_OverlayTopologyError("Found null edge in ring"))
        edge = ne
        #-- `ids` stays OPEN: the final dest is the start origin, already pushed
        edge == ring.start_edge && (_ring_add!(pts, node_point(ctx.arr, dest)); break)
        push!(ids, Int32(dest))
        _ring_add!(pts, node_point(ctx.arr, dest))
    end
    #-- the last dest is the start origin, so pts is already closed; be defensive
    pts[end] == pts[1] || push!(pts, pts[1])

    ring.ring_pts = pts
    ring.node_ids = ids
    ring.is_hole = _ring_is_ccw_exact(ctx, ids)
    ring.bbox = _ring_bbox(pts)
    return nothing
end

#=
## Shell-vs-hole from exact kernel predicates (design §0)

Result assembly is a decision, so it is decided on the arrangement, not on its
emitted image. `_ring_is_ccw_exact` answers "is this minimal ring wound CCW?"
from the ring's NODE IDS — the symbolic, exact identities — never from
`ring_pts`, which is `node_point` output: rounded once at emission (design
§2.6). The distinction is not academic. A ring whose exact width is a fraction
of an ULP rounds to a self-touching or inverted image, and `Orientation.isCCW`
over that image reports the opposite role; the ring is then filed as a second
shell or as a free hole no shell contains, which is where a family of
`OverlayTopologyError`s came from.

Deciding the ROLE exactly does not close that family, and the comment above used
to imply it did. A sub-ULP ring now gets its role right and is then emitted, and
its emitted image is an invalid ring — the failure moves from a throw to invalid
output, which is the same contract violation wearing different clothes. Closing
it needs the ring dropped rather than merely classified, which is what
`_ring_is_subgrid` does; see the emission-grid discussion above. (Sub-ULP rings
are also not the only source of `OverlayTopologyError` — the labeller has an
independent one, tracked separately.)

Both manifolds take the sign of the ring's exact signed area — the planar
shoelace, the spherical geodesic curvature. Neither reads an output coordinate.
=#

#=
Planar: the sign of the doubled signed area over the exact node coordinates,
as a certified Float64 filter with an exact `Rational{BigInt}` fallback — the
same filter/escalate shape every other predicate in the port uses.

The filter runs on `node_point`, which is the CORRECTLY ROUNDED exact node
(design §2.6: certified double-double, else the rounded rational), so each
coordinate carries at most ½ ulp of error. The bound below covers both that
perturbation (`cmax * len`, a coordinate displacement acting on the ring's
`L1` extent) and the accumulated summation error (`n * mag`). Clearing it
certifies the sign of the EXACT area, not merely of the emitted one; failing
it escalates.

The signed area is used rather than JTS's extreme-vertex "cap" test because
its sign is the ring's net winding whether or not the ring is simple — the cap
test presumes simplicity, which holds for arrangement minimal rings but is one
more premise to carry. Translating to the first vertex costs nothing and keeps
the filter well-conditioned on the corpus's ~1e6-magnitude coordinates.
=#
function _ring_is_ccw_exact(ctx::_PolyBuilderCtx{<:Planar}, ids::Vector{Int32})
    length(ids) < 3 && return false
    (area, certain) = _ring_area2_filter(ctx.arr, ids)
    certain && return area > 0
    return _ring_area2_exact(ctx.arr, ids) > 0
end

# The accumulated Float64 doubled signed area and its error bound. Both the sign
# test above and the magnitude test in `_ring_width_below` are read off these.
function _ring_area2_bounded(arr, ids::Vector{Int32})
    n = length(ids)
    ox, oy = node_point(arr, ids[1])
    acc = 0.0; mag = 0.0; len = 0.0; cmax = max(abs(ox), abs(oy))
    ax = 0.0; ay = 0.0                       # translated previous vertex
    @inbounds for i in 2:n
        (px, py) = node_point(arr, ids[i])
        cmax = max(cmax, abs(px), abs(py))
        bx = px - ox; by = py - oy
        t1 = ax * by; t2 = ay * bx
        acc += t1 - t2
        mag += abs(t1) + abs(t2)
        len += abs(bx) + abs(by)
        ax = bx; ay = by
    end
    #-- the closing term (last → first) has both translated factors zero
    u = 0.5 * eps(Float64)
    bound = 8 * u * (n * mag + cmax * len)
    return (acc, bound)
end

# Exact doubled signed area: the plain shoelace over the arrangement's exact
# node coordinates (`_exact_node_point` — the input vertex for a vertex node,
# `_exact_crossing_point` for a crossing). Rational arithmetic on Float64-derived
# values is exact, so the returned sign is the ring's true winding.
function _ring_area2_exact(arr, ids::Vector{Int32})
    n = length(ids)
    acc = zero(Rational{BigInt})
    (ax, ay) = _exact_node_point(arr.nodes.keys[ids[1]])
    (x1, y1) = (ax, ay)
    @inbounds for i in 2:n
        (bx, by) = _exact_node_point(arr.nodes.keys[ids[i]])
        acc += ax * by - ay * bx
        ax = bx; ay = by
    end
    return acc + (ax * y1 - ay * x1)
end

#=
Spherical: the loop's geodesic curvature (`_ring_is_ccw(::Spherical, …)`, our
port of S2 `GetCurvature`) over the exact node directions. The turning-angle
formulation is apex-free — each term sees only three adjacent vertices, and its
sign comes from the exact `Sign` predicate — so it is the formulation this
codebase has already settled on for spherical orientation, and the planar
extreme-vertex cap is deliberately NOT ported to it.

What changes here is only its input: each vertex is the node's own direction —
the ingested unit vector for a vertex node, the exact crossing direction
`±(na × nb)` for a crossing node — instead of the emitted `(lon, lat)` converted
back through trigonometry. Distinct nodes therefore stay distinct, which
`_prune_loop_degeneracies` would otherwise collapse.
=#
_ring_is_ccw_exact(ctx::_PolyBuilderCtx{<:Spherical}, ids::Vector{Int32}) =
    _ring_is_ccw(ctx.m, [_node_kernel_point(ctx, id) for id in ids]; exact = ctx.exact)

# The exact sphere position of node `id` as a unit kernel point, memoized per
# builder (a node is shared by every ring through it, and the exact crossing
# direction is `Rational{BigInt}` work).
function _node_kernel_point(ctx::_PolyBuilderCtx{<:Spherical, P}, id::Integer) where {P}
    if isempty(ctx.kernel_ok)
        n = num_nodes(ctx.arr)
        ctx.kernel_cache = Vector{P}(undef, n)
        ctx.kernel_ok = fill(false, n)
    end
    i = Int(id)
    @inbounds ctx.kernel_ok[i] && return ctx.kernel_cache[i]
    k = ctx.arr.nodes.keys[i]
    p = if k.is_crossing
        d = _sph_crossing_dir(booltype(ctx.exact), k)
        rk_normalize_usp(UnitSphericalPoint(Float64(d[1]), Float64(d[2]), Float64(d[3])))
    else
        k.pt
    end
    @inbounds ctx.kernel_cache[i] = p
    @inbounds ctx.kernel_ok[i] = true
    return p
end

#=
Spherical counterpart of `_ring_width_below`, same one-sided contract.

The area used is the CHORDAL one — the area of the flat polygon spanned by the
ring's kernel directions in R³, `½|Σ (vᵢ−v₁) × (vᵢ₊₁−v₁)|` — rather than the
spherical area, for one reason: it is a polynomial in the coordinates, so it
escalates to exact `Rational{BigInt}` the same way the planar shoelace does. The
spherical area is transcendental (Girard, or the Van Oosterom–Strackee
`2·atan((a×b)·c, 1+a·b+b·c+c·a)`) and has no exact rational form to escalate to;
it would need extended-precision floating point and a bespoke error analysis.

The substitution costs nothing here. Chord and arc agree to relative `d²/24` for
a separation `d`, and this test only ever fires on rings of mean width below
`4u ≈ 5e-16` radians, whose chordal diameter `d` in the corpus is ~1e-5 — a
relative discrepancy of ~4e-12. It cannot fire spuriously on a large ring either:
for any ring spanning an appreciable part of the sphere both the chordal area and
the extent below are O(1), so the ratio is O(1) and sits sixteen orders of
magnitude above `4u`.

`P` is the chordal perimeter and enters exactly as it does on the plane: a
Float64 filter with an error bound, escalating to a certified rational lower
bound when that cannot decide. The kernel points ARE Float64 triples, so the
escalated perimeter carries no endpoint displacement — only the square roots,
which `_sqrt_lo` backs off.
=#
function _ring_width_below(ctx::_PolyBuilderCtx{<:Spherical}, ids::Vector{Int32}, thr::Float64)
    n = length(ids)
    ks = [_node_kernel_point(ctx, id) for id in ids]
    #-- chordal perimeter, Float64, with its rounding + endpoint error bound
    per = 0.0
    @inbounds for i in 1:n
        a = ks[i]; b = ks[i == n ? 1 : i + 1]
        f1 = b[1] - a[1]; f2 = b[2] - a[2]; f3 = b[3] - a[3]
        per += sqrt(f1 * f1 + f2 * f2 + f3 * f3)
    end
    per_err = 4 * eps(Float64) * n * (1.0 + per)     # unit vectors: |coord| ≤ 1
    #-- Float64 filter on |S| where S = Σ (vᵢ−v₁) × (vᵢ₊₁−v₁), so 2·A_chord = |S|
    o = ks[1]
    s1 = 0.0; s2 = 0.0; s3 = 0.0
    m1 = 0.0; m2 = 0.0; m3 = 0.0                     # Σ|individual products|
    L = 0.0                                          # max translated |component|
    ax = 0.0; ay = 0.0; az = 0.0
    @inbounds for i in 2:n
        p = ks[i]
        bx = p[1] - o[1]; by = p[2] - o[2]; bz = p[3] - o[3]
        L = max(L, abs(bx), abs(by), abs(bz))
        t1 = ay * bz; t2 = az * by
        t3 = az * bx; t4 = ax * bz
        t5 = ax * by; t6 = ay * bx
        s1 += t1 - t2; s2 += t3 - t4; s3 += t5 - t6
        m1 += abs(t1) + abs(t2); m2 += abs(t3) + abs(t4); m3 += abs(t5) + abs(t6)
        ax = bx; ay = by; az = bz
    end
    #-- the closing term (last → first) has both translated factors zero
    e = eps(Float64)
    b1 = 8 * e * (n * m1 + n * L * L)
    b2 = 8 * e * (n * m2 + n * L * L)
    b3 = 8 * e * (n * m3 + n * L * L)
    sn = sqrt(s1 * s1 + s2 * s2 + s3 * s3)
    bn = b1 + b2 + b3                                # ≥ ‖(b1,b2,b3)‖, so ≥ |Δ‖S‖|
    if per > per_err
        sn - bn > thr * (per + per_err) && return false
        sn + bn < thr * (per - per_err) && return true
    end
    #-- escalate: S exactly, P as a certified rational lower bound
    R = Rational{BigInt}
    o1 = R(o[1]); o2 = R(o[2]); o3 = R(o[3])
    q1 = zero(R); q2 = zero(R); q3 = zero(R)
    rax = zero(R); ray = zero(R); raz = zero(R)
    @inbounds for i in 2:n
        p = ks[i]
        rbx = R(p[1]) - o1; rby = R(p[2]) - o2; rbz = R(p[3]) - o3
        q1 += ray * rbz - raz * rby
        q2 += raz * rbx - rax * rbz
        q3 += rax * rby - ray * rbx
        rax = rbx; ray = rby; raz = rbz
    end
    plo = zero(R)
    @inbounds for i in 1:n
        a = ks[i]; b = ks[i == n ? 1 : i + 1]
        f1 = R(b[1]) - R(a[1]); f2 = R(b[2]) - R(a[2]); f3 = R(b[3]) - R(a[3])
        plo += _sqrt_lo(f1 * f1 + f2 * f2 + f3 * f3)
    end
    #-- |S| < thr·P  <=>  |S|² < (thr·P)², both sides non-negative
    rhs = R(thr) * plo
    return q1 * q1 + q2 * q2 + q3 * q3 < rhs * rhs
end

function _ring_bbox(pts::Vector{Tuple{Float64, Float64}})
    xmin = xmax = pts[1][1]
    ymin = ymax = pts[1][2]
    for p in pts
        xmin = min(xmin, p[1]); xmax = max(xmax, p[1])
        ymin = min(ymin, p[2]); ymax = max(ymax, p[2])
    end
    return (xmin, xmax, ymin, ymax)
end

# ## Hole containment (ports of `OverlayEdgeRing.locate` / `contains` / …)
#
# KNOWN GAP, tracked separately from the shell/hole role above: this half of
# result assembly still decides on `ring_pts`, i.e. on emitted coordinates. The
# ray crossing itself is robust (`rk_orient`, never naive even-odd), but its
# inputs are rounded, and so are the bbox/RTree prefilters that reject candidate
# shells. Making it exact needs a point-in-ring predicate over the arrangement's
# node coordinates on both manifolds (`Rational{BigInt}` crossing counts on the
# plane; a rational-direction crossing-parity scan on the sphere), plus outward
# padding of the two float prefilters so they cannot reject a genuinely contained
# hole, plus a way to index rational coordinates (or to demote the existing
# Float64 y-interval index to a candidate-edge filter). Measured over the robust
# corpus, all 405 free-hole placements are decided by a strictly-INTERIOR hole
# vertex — none falls through the all-BOUNDARY path — so nothing currently turns
# on it.

# Lazily builds an indexed point-in-area locator over this ring and locates `p`
# (design §3 amendment 5: robust ray crossing via `rk_orient`, over the ring's
# emitted coordinates — never naive even-odd).
function _ring_locate(ctx, ring::_OverlayEdgeRing, p)
    if ring.locator === nothing
        ring.locator = IndexedPointInAreaLocator(ctx.m, GI.Polygon([ring.ring_pts]);
                                                 exact = ctx.exact)
    end
    return locate(ring.locator, p)
end

# Whether `shell` contains `hole` (port of `contains` + `isPointInOrOut`). On the
# plane a bounding-box reject prefilters before the point tests; on the sphere the
# lon/lat box is unreliable near the poles/antimeridian, so containment is decided
# purely by the point tests (free holes are rare, this path is cold).
#
# Adaptation to the non-self-noding substrate: JTS uses `containsProperly` here,
# because a hole touching its own shell would have been connected into one maximal
# ring (JTS nodes A against itself). This substrate does NOT self-node a single
# input (design §2.2), so such a hole surfaces as a disconnected free hole whose
# bbox touches its shell's; the prefilter must therefore be non-strict
# (`_bbox_contains`), letting the point-in-area test — the real decision — run.
_ring_contains(ctx::_PolyBuilderCtx{<:Planar}, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing) =
    _bbox_contains(shell.bbox, hole.bbox) && _is_point_in_or_out(ctx, shell, hole)
_ring_contains(ctx::_PolyBuilderCtx{<:Spherical}, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing) =
    _is_point_in_or_out(ctx, shell, hole)

function _is_point_in_or_out(ctx, shell::_OverlayEdgeRing, hole::_OverlayEdgeRing)
    for p in hole.ring_pts
        loc = _ring_locate(ctx, shell, p)
        loc == LOC_INTERIOR && return true
        loc == LOC_EXTERIOR && return false
        #-- LOC_BOUNDARY: inconclusive, keep checking
    end
    return false
end

@inline _bbox_contains(outer::NTuple{4, Float64}, inner::NTuple{4, Float64}) =
    outer[1] <= inner[1] && outer[2] >= inner[2] && outer[3] <= inner[3] && outer[4] >= inner[4]

# Port of `findEdgeRingContaining`: the innermost (smallest-envelope) shell in
# `candidates` that contains `hole`, or `0`.
function _find_edge_ring_containing(ctx, hole::_OverlayEdgeRing, candidates)
    min_containing = Int32(0)
    for sh in candidates
        shell = ctx.edge_rings[sh]
        if _ring_contains(ctx, shell, hole)
            if min_containing == 0 ||
               _bbox_contains(ctx.edge_rings[min_containing].bbox, shell.bbox)
                min_containing = Int32(sh)
            end
        end
    end
    return min_containing
end
