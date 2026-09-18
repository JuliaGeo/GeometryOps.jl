# NOTE: This functionality is experimental and may change at any time.

# # Monotone chains as noding index items (design §2.3)
#
#=
Port of JTS `MCIndexNoder` / `index.chain.MonotoneChain`: instead of indexing
one extent per segment, index one extent per *monotone chain* — a maximal run of
consecutive segments whose direction stays in one quadrant. A chain is weakly
monotone in both coordinates, so its envelope, and the envelope of any
contiguous section of it, is the box of that section's two endpoints. Candidate
segment pairs then come from a mutual binary descent over two chains
(`_chain_segment_pairs`), which rejects a whole section pair on one endpoint-box
test at every level.

Two properties do the work on offset-curve linework:

  * **Far fewer index items.** A fillet quarter-arc is one chain whatever its
    segment count; a smooth 10 000-vertex circle is four chains. The R-tree's
    intra-leaf all-pairs work is quadratic in items per leaf, so this is where
    the traversal time goes.
  * **Within-chain pairs are never enumerated.** Two segments of one monotone
    chain cannot record anything: both coordinates are monotone along the chain,
    so sections that are not adjacent have disjoint interiors, and adjacent ones
    meet only at the vertex between them — an endpoint of both, which
    `_classify_pair!`'s guards drop. That is the same reason JTS never compares a
    chain with itself.

Chains are planar-only: `_segment_extent(::Spherical, p, q)` is a great-circle
arc extent, which is the box of the arc *between* two points, not of the polyline
joining them, so it does not bound a chain.
=#

# JTS `Quadrant.quadrant`: NE=0, NW=1, SW=2, SE=3, axis directions belonging to
# the `dx >= 0` / `dy >= 0` side. A zero-length direction lands in NE and is
# never consulted (`_chain_end` skips zero-length segments, as JTS does).
@inline function _seg_quadrant(p, q)
    dx = GI.x(q) - GI.x(p)
    dy = GI.y(q) - GI.y(p)
    return dx >= 0 ? (dy >= 0 ? 0 : 3) : (dy >= 0 ? 1 : 2)
end

# Port of `MonotoneChainBuilder.findChainEnd`: the index of the last point of the
# chain starting at point `start`. Zero-length segments join the chain but never
# set or break its quadrant.
function _chain_end(pts, start::Int)
    n = length(pts)
    safe = start
    while safe < n && pts[safe] == pts[safe + 1]
        safe += 1
    end
    safe >= n && return n
    q = _seg_quadrant(pts[safe], pts[safe + 1])
    last = start + 1
    while last < n
        if pts[last] != pts[last + 1] && _seg_quadrant(pts[last], pts[last + 1]) != q
            break
        end
        last += 1
    end
    return last
end

#=
The chain counterpart of `_relate_edge_index`: an `Unsorted` R-tree over chain
envelopes whose leaf data is `(string index, first point, last point)`. Wrapped
in `ChainIndex` so `_collect_self_pairs!` can dispatch on it — a segment index
and a chain index are the same `RTree` type and differ only in what a leaf means.
=#
struct ChainIndex{TR}
    tree::TR
    nchains::Int
    nsegments::Int
end

function _chain_extent_table(m::Planar, ss_list)
    extents = _segment_extent_type(m)[]
    owners = NTuple{3, Int}[]
    nseg = _total_segment_count(ss_list)
    sizehint!(extents, nseg ÷ 4 + 8)
    sizehint!(owners, nseg ÷ 4 + 8)
    for (si, ss) in enumerate(ss_list)
        _push_chains!(m, extents, owners, si, ss.pts)
    end
    return extents, owners, nseg
end

# Function barrier: one dispatch per string, as in `_push_segment_extents!`.
function _push_chains!(m::Planar, extents::Vector, owners::Vector, si::Int, pts::Vector)
    n = length(pts)
    lo = 1
    while lo < n
        hi = _chain_end(pts, lo)
        push!(owners, (si, lo, hi))
        push!(extents, _segment_extent(m, pts[lo], pts[hi]))
        lo = hi
    end
    return nothing
end

function _relate_chain_index(m::Planar, ss_list)
    extents, owners, nseg = _chain_extent_table(m, ss_list)
    isempty(extents) && return nothing
    return ChainIndex(RTree(Unsorted(), owners; extents, nodecapacity = 16),
                      length(owners), nseg)
end

#=
The index the self-noding pass traverses: one leaf per monotone chain, always.

Chains buy exactly one thing — the `nsegments - nchains` within-chain adjacent
pairs are never enumerated — and charge for two: a chain envelope is looser than
a segment's, and each surviving pair arrives through a descent rather than
straight off the index. The buy wins whenever the curve turns outward more often
than it reverses, since a run that never reverses is one chain however long. An
offset curve that reverses at every vertex is the shape to watch: its chains are
one segment each and the descent is then pure overhead.
=#
function _noding_index(m::Planar, ss_list)
    _total_segment_count(ss_list) == 0 && return nothing
    return _relate_chain_index(m, ss_list)
end

_noding_index(m::Manifold, ss_list) = _relate_edge_index(m, ss_list)

#=
Port of `MonotoneChain.computeOverlaps`: mutual binary descent over two chain
sections, `f(k0, k1)` once per segment pair whose own boxes overlap. `s`/`e` are
point indices, so a section spans segments `s:(e - 1)`. The caller has already
tested the two whole-chain envelopes — they are the extents the index holds.

Two departures from the Java, both about not paying for the descent where it
cannot pay for itself. Offset-curve linework is full of two-segment chains (an
inside turn reverses direction, which ends a chain), so the common case is a
section pair small enough that halving it costs more than testing its segment
pairs outright:

  * a section pair of at most 2×2 segments is enumerated directly, one
    endpoint-box test per segment pair — the same test, and the same count, the
    per-segment index would have applied;
  * a 1×1 section is handed over with no test at all, because the test that
    admitted it was on exactly that box.

The box test therefore runs at every level including the leaf, where JTS stops
testing and hands the pair over unconditionally. That makes a pair reach
`_classify_pair!` exactly when its two segment boxes overlap — the same set the
per-segment index delivers, so the two enumerations are directly comparable and
the classifier sees identical input.
=#
@inline function _chain_segment_pairs(f::F, pa, s0::Int, e0::Int, pb, s1::Int, e1::Int) where {F}
    n0 = e0 - s0
    n1 = e1 - s1
    if (n0 == 1) & (n1 == 1)
        f(s0, s1)
    elseif (n0 <= 2) & (n1 <= 2)
        @inbounds for k0 in s0:(e0 - 1), k1 in s1:(e1 - 1)
            _endpoint_boxes_overlap(pa[k0], pa[k0 + 1], pb[k1], pb[k1 + 1]) && f(k0, k1)
        end
    else
        _chain_descend(f, pa, s0, e0, pb, s1, e1)
    end
    return nothing
end

# The recursive half, kept out of line so the short-chain cases above inline into
# the index traversal's leaf callback.
function _chain_descend(f::F, pa, s0::Int, e0::Int, pb, s1::Int, e1::Int) where {F}
    m0 = (s0 + e0) >> 1
    m1 = (s1 + e1) >> 1
    if s0 < m0
        s1 < m1 && _chain_try(f, pa, s0, m0, pb, s1, m1)
        m1 < e1 && _chain_try(f, pa, s0, m0, pb, m1, e1)
    end
    if m0 < e0
        s1 < m1 && _chain_try(f, pa, m0, e0, pb, s1, m1)
        m1 < e1 && _chain_try(f, pa, m0, e0, pb, m1, e1)
    end
    return nothing
end

@inline function _chain_try(f::F, pa, s0::Int, e0::Int, pb, s1::Int, e1::Int) where {F}
    @inbounds _endpoint_boxes_overlap(pa[s0], pa[e0], pb[s1], pb[e1]) &&
        _chain_segment_pairs(f, pa, s0, e0, pb, s1, e1)
    return nothing
end

# Closed-interval overlap of the boxes of (p0, p1) and (q0, q1) — the same test
# `Extents.intersects` applies to the stored extents.
@inline function _endpoint_boxes_overlap(p0, p1, q0, q1)
    plo, phi = minmax(GI.x(p0), GI.x(p1))
    qlo, qhi = minmax(GI.x(q0), GI.x(q1))
    ((plo > qhi) | (phi < qlo)) && return false
    plo, phi = minmax(GI.y(p0), GI.y(p1))
    qlo, qhi = minmax(GI.y(q0), GI.y(q1))
    return !((plo > qhi) | (phi < qlo))
end
