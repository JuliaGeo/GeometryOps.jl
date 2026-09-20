# NOTE: This functionality is experimental and may change at any time.

# # Stage 1 — collect (design §2.3)
#
# Enumerate candidate segment pairs (A × B) through the reused RelateNG edge
# index and classify each with the exact kernel (`rk_classify_intersection`),
# recording symbolic node keys onto the parent segments. No intersection
# coordinate is ever constructed here.
#
# `seg_nodes` accumulates, per parent segment `(string_idx, seg_idx)`, the ids of
# nodes lying strictly in that segment's interior. `string_idx` is global into
# the arrangement's `segstrings`: A strings occupy `1:na`, B strings `na+1:end`.

#=
`seg_nodes` is a flat, push-only `(string, segment, node)` list, sorted once
after collect (`_merge_coincident_nodes!`) and then read by `split.jl` in one
lockstep walk with its own segment loop. The keyed-`Dict` shape it replaces paid
a hash on every record *and* a hash on every segment of every string at split
time, including the overwhelming majority that carry no interior node at all.
=#
@inline function _record_interior!(seg_nodes::Vector{NTuple{3, Int32}},
        string_idx::Int32, seg_idx::Int32, nid::Int32)
    push!(seg_nodes, (string_idx, seg_idx, nid))
    return nothing
end

function _collect_crossings!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ssa::AbstractVector{RelateSegmentString{P}},
        ssb::AbstractVector{RelateSegmentString{P}}, na::Int32;
        exact = True(), tree_a = nothing, tree_b = nothing,
        clip_a = nothing, clip_b = nothing, self_node_all::Bool = false) where {P}
    ta = tree_a === nothing ? _relate_edge_index(m, ssa) : tree_a
    tb = tree_b === nothing ? _relate_edge_index(m, ssb) : tree_b
    #-- a chain index (`chains.jl`) indexes runs of segments, not segments, so it
    #-- only ever reaches the self-noding pass below
    (ta isa ChainIndex || tb isa ChainIndex) && !isempty(ssb) &&
        throw(ArgumentError("a monotone-chain index cannot drive the A x B noding pass"))
    if !(ta === nothing || tb === nothing)
        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, ta, tb) do ia, ib
            (sa, ka) = ta.data[ia]
            (sb, kb) = tb.data[ib]
            _classify_pair!(m, table, seg_nodes, ssa, sa, Int32(sa), Int32(ka),
                            ssb, sb, na + Int32(sb), Int32(kb); exact)
            return nothing
        end
    end
    #-- design §2.2 amendment: each input must additionally be self-noded, in two
    #-- passes of different scope (see below). Both are restricted to the side's
    #-- clip box when it has one — see "Self-noding under clip pruning".
    _collect_self_crossings!(m, table, seg_nodes, ssa, Int32(0);
                             exact, clip = clip_a, all_strings = self_node_all, tree = ta)
    _collect_self_crossings!(m, table, seg_nodes, ssb, na;
                             exact, clip = clip_b, all_strings = self_node_all, tree = tb)
    #-- `self_node_all` has already run the all-pairs pass over every string, which
    #-- subsumes the vertex pass (it classifies each segment pair exactly, vertex
    #-- incidences included) — the same reason the vertex pass skips an all-linear side
    if !self_node_all
        _collect_self_vertex_nodes!(m, table, seg_nodes, ssa, Int32(0), ta; exact, clip = clip_a)
        _collect_self_vertex_nodes!(m, table, seg_nodes, ssb, na, tb; exact, clip = clip_b)
    end
    return nothing
end

#=
## Self-noding under clip pruning

(Scope, as in split.jl: a clip box exists only for a PLANAR intersection or
difference. `clip` is `nothing` on the sphere and for every union, and the
filters below then reduce to their `::Nothing` methods.)

The A×B pass above needs no clip filter: a pruned A segment's bbox misses
`env(B)` (that is what pruning to `env(A) ∩ env(B)`, or to `env(A)`, means), so
it never shares an extent with any B segment and the dual traversal already
never enumerates it. Self-noding is different — a side's own segments overlap
each other constantly, there is no pruning to be had from the data, and the pass
is O(all segments) whatever the box. Measured on Natural Earth 10 m Brazil
against a 2° box, it was 947 µs of a 1.83 ms query with every other stage already
pruned.

So both self-noding passes take the side's clip box and restrict the traversal
to it: a tree node whose extent misses the box is never descended into, and a
segment whose extent misses the box is never paired. That is exactly the same
kept/pruned test `split.jl` applies (`_segment_extent(::Planar, p, q)` IS the
endpoint bbox), so the two stages agree segment for segment.

**What this drops, and why it is sound.** Only pairs with at least one PRUNED
member. A pruned segment lies wholly outside the box, so any node such a pair
would create lies outside the box too. Three cases:

  * both pruned — neither emits an edge, so the node would reference nothing;
  * a pruned segment's vertex inside a KEPT segment — the kept segment is then
    left unsplit there;
  * a kept segment's vertex inside a PRUNED segment — the pruned segment emits
    no edge either way, so only the node's existence is lost.

The middle case is the one that changes the graph, and it cannot change the
result. The A×B pass is untouched, so every kept segment is still split at every
incidence with the OTHER input; between consecutive nodes a kept edge therefore
has a uniform location with respect to the other input, and its label — from the
ring's own `depth_delta` plus that uniform location — is correct along its whole
length. An edge marked in-result consequently lies wholly on the result boundary,
hence inside the result's closure, hence inside the box; an edge spanning a
dropped same-side node runs through a point outside the box and so is never
marked. What is lost is only the subdivision of NON-result linework at
same-side touch points outside the box.

Two invariants survive verbatim. **Every node inside the box keeps its full
star**, and the arrangement restricted to the box is identical to the unpruned
one: any incidence at a point in the box involves only segments whose bbox meets
the box, and those are all still in the traversal. And coincident segments have
identical extents, so they prune together and are noded in lockstep — the edge
merger never sees a half-split pair. Ring assembly, which only ever links
in-result edges, is therefore reading exactly the graph it would have read.

The stars that DO thin are handled by `arr.truncated`, which `split.jl` derives
positionally (a vertex node outside the box) rather than from a record of what
each stage dropped — which is what lets this pass drop nodes without having to
report them.
=#

#=
## Self-noding (design §2.2 amendment)

The arrangement invariant is that no node lies strictly inside a noded edge
(design §2.1). The A×B pass alone does not establish it: a node lying in the
interior of a segment of the *same* input as one of the pair is never created by
an A×B classification, so nothing splits that segment. `TestOverlayLA.xml` case 2
is the canonical instance — one MULTILINESTRING component runs along a polygon
hole boundary while a second component of the same MULTILINESTRING crosses it,
so the hole edge is split at the two crossings but the collinear line component
is not, the two no longer share a node pair, and the edge merger never pairs
them.

Self-noding therefore runs per side, in two passes of *different scope*, because
what a valid input can do to itself depends on its dimension:

- LINEAR strings get the full all-pairs pass below. A valid MultiLineString
  carries no self-noding guarantee whatsoever: components may properly cross,
  overlap collinearly, and retrace each other.
- Every string's VERTICES are additionally noded against the side's own segments
  (`_collect_self_vertex_nodes!`). For a valid AREAL input that is the whole of
  what self-noding can find, so the expensive all-pairs pass is not run for it.

The areal argument. Rings of a valid polygonal geometry never cross and never
overlap: any two of them meet in at most finitely many points, so a segment pair
drawn from one input intersects either in nothing or in a single point. A single
point shared by two segments that do not properly cross is an endpoint of at
least one of them — i.e. an input vertex. So the only self-incidence a valid area
can carry is *vertex of one ring inside the interior of another ring's segment*
(a hole apex on its shell, or two multipolygon components meeting), which is
exactly what the vertex pass records. Note this is not a claim that areas need no
self-noding: they do, and the missing node is a real defect — the fuzz suite's
`hole apex on the shell edge (self-touching input)` class, 4 broken tests before
this pass existed.

That argument was checked, not assumed. Over 5 445 input pairs — Natural Earth
110 m and 10 m country neighbours (planar and spherical) and shifted-self copies,
600 cases from each of the eight fuzz generators, and every case of the JTS
`overlay` and `overlay_robust` XML corpora — the per-segment interior node sets
produced by this vertex pass and by the full all-pairs pass run on every string
are IDENTICAL on 5 444. The one input they differ on, `TestOverlayMisc.xml`
case 5 (GEOS ticket 737), is not valid: `isValidReason` reports a
self-intersection and its rings carry four proper self-crossings, which is
outside the engine's contract. Both scopes throw `OverlayTopologyError` on it,
before and after this change.

Cost, over 29 Natural Earth 10 m country pairs (16k-53k points per pair, 884k
total), arrangement build time, best of five:

    0.220 s   no areal self-noding (what this branch did before)
    0.346 s   the targeted pass here
    0.474 s   the full all-pairs pass on every string, `dual(t, t)` enumeration

Read the breakdown before attributing that to the scope. Timing the self-noding
pass alone over the 46 sides of that corpus:

    0.258 s   `dual(t, t)` enumeration + full classification
    0.137 s   `_self_pair_search` enumeration + full classification
    0.131 s   `_self_pair_search` enumeration + vertex-only test  (ships)

So nearly all of the saving is the unordered-pair enumeration, which the full
pass could have had too; narrowing the scope to vertices is worth ~5% on top.
The scope is still the right one — it is provably sufficient on valid areal
input, it is what the argument above licenses, and it never fabricates a
crossing node on an input that should not have one — but it is not where the
time went.
=#
#=
`all_strings` widens the scope from LINEAR strings to every string on the side,
and is what a caller asks for when its input is area-contributing but *not* a
valid area — the N-ary winding overlay's self-intersecting offset curves. The
dimension test above is a soundness argument about VALID input; linework that
carries no validity guarantee needs the all-pairs pass whatever its dimension,
and the caller is the only party that knows which it has.

With every string selected, `sub` is the caller's own list, so `tree` — the index
the A x B pass already built over exactly that list — is reused instead of a
second one being built over a copy of it.
=#
function _collect_self_crossings!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ss::AbstractVector{RelateSegmentString{P}}, off::Int32;
        exact, clip = nothing, all_strings::Bool = false, tree = nothing) where {P}
    if all_strings
        #-- fewer than two segments in total: no pair to classify
        sum(s -> length(s.pts) - 1, ss; init = 0) < 2 && return nothing
        t = tree === nothing ? _noding_index(m, ss) : tree
        t === nothing && return nothing
        return _collect_self_pairs!(m, table, seg_nodes, ss, off, nothing, t, clip; exact)
    end
    lin = Int32[Int32(i) for i in eachindex(ss) if ss[i].dim == DIM_L]
    isempty(lin) && return nothing
    sub = [ss[i] for i in lin]
    sum(s -> length(s.pts) - 1, sub; init = 0) < 2 && return nothing
    t = _relate_edge_index(m, sub)
    t === nothing && return nothing
    return _collect_self_pairs!(m, table, seg_nodes, sub, off, lin, t, clip; exact)
end

# `local_to_global[s]` maps a position in `ss` to its arrangement-global string
# index; `nothing` when `ss` IS the side's own list and the map is the identity.
@inline _self_string_idx(::Nothing, s::Integer) = Int32(s)
@inline _self_string_idx(map::Vector{Int32}, s::Integer) = @inbounds map[s]

#=
Chain-indexed enumeration (`chains.jl`): the tree's leaves are monotone chains,
so a leaf pair is a *chain* pair and the segment pairs inside it come from the
mutual binary descent. Clip pruning is not offered here — the chain envelope
keeps segments its members would individually have been pruned by, which would
change which same-side nodes are dropped; the only caller (`winding_overlay.jl`)
never clips.
=#
function _collect_self_pairs!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ss::AbstractVector{RelateSegmentString{P}}, off::Int32,
        local_to_global, ci::ChainIndex, clip; exact) where {P}
    clip === nothing ||
        throw(ArgumentError("chain-indexed self-noding does not support clip pruning"))
    t = ci.tree
    _self_pair_search(m, t, nothing) do i1, i2
        (s1, lo1, hi1) = t.data[i1]
        (s2, lo2, hi2) = t.data[i2]
        p1 = ss[s1].pts
        p2 = ss[s2].pts
        g1 = off + _self_string_idx(local_to_global, s1)
        g2 = off + _self_string_idx(local_to_global, s2)
        _chain_segment_pairs(p1, lo1, hi1, p2, lo2, hi2) do k1, k2
            _classify_pair!(m, table, seg_nodes, ss, s1, g1, Int32(k1),
                            ss, s2, g2, Int32(k2); exact)
            return nothing
        end
        return nothing
    end
    return nothing
end

function _collect_self_pairs!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ss::AbstractVector{RelateSegmentString{P}}, off::Int32,
        local_to_global, t, clip; exact) where {P}
    _self_pair_search(m, t, clip) do i1, i2
        (s1, k1) = t.data[i1]
        (s2, k2) = t.data[i2]
        _classify_pair!(m, table, seg_nodes,
                        ss, s1, off + _self_string_idx(local_to_global, s1), Int32(k1),
                        ss, s2, off + _self_string_idx(local_to_global, s2), Int32(k2); exact)
        return nothing
    end
    return nothing
end

#=
The targeted half of self-noding: one side's vertices against that side's own
segments, recording a vertex as an interior node of every segment whose
*interior* it lies in, and recording nothing else. Sufficient scope for valid
areal input — argued and measured above.

Enumeration reuses the side's existing segment index — the one the A×B pass
already built, or the one the caller supplied — so no second index is
constructed, and it walks that index against itself with `_self_pair_search`,
visiting each unordered segment pair once. A vertex in a segment's interior is
an endpoint of some segment whose extent therefore overlaps that segment's, so
the pair is enumerated and the incidence is found; testing the four
vertex/segment combinations at each pair is what makes this a vertex pass rather
than an intersection pass. `tree.data[i]` is `(string index within `ss`, segment
index)`, hence the `off` to reach the arrangement-global string index.

The membership test is `rk_point_on_segment`, the exact kernel predicate, minus
the two endpoints — no tolerance, no distance, in keeping with design §0. A
vertex equal to a segment endpoint is rejected before the predicate runs: it is
already a node of that segment and needs no split, and rejecting it early keeps
every adjacent-segment and shared-ring-vertex pair off the predicate entirely.
=#
function _collect_self_vertex_nodes!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ss::AbstractVector{RelateSegmentString{P}}, off::Int32, tree;
        exact, clip = nothing) where {P}
    tree === nothing && return nothing
    #-- an all-linear side is fully covered by the pass above
    all(s -> s.dim == DIM_L, ss) && return nothing
    _self_pair_search(m, tree, clip) do i1, i2
        (s1, k1) = tree.data[i1]
        (s2, k2) = tree.data[i2]
        p1 = ss[s1].pts; p2 = ss[s2].pts
        a0 = p1[k1]; a1 = p1[k1 + 1]
        b0 = p2[k2]; b1 = p2[k2 + 1]
        _record_vertex_on_segment!(m, table, seg_nodes, a0, b0, b1, off + Int32(s2), Int32(k2); exact)
        _record_vertex_on_segment!(m, table, seg_nodes, a1, b0, b1, off + Int32(s2), Int32(k2); exact)
        _record_vertex_on_segment!(m, table, seg_nodes, b0, a0, a1, off + Int32(s1), Int32(k1); exact)
        _record_vertex_on_segment!(m, table, seg_nodes, b1, a0, a1, off + Int32(s1), Int32(k1); exact)
        return nothing
    end
    return nothing
end

# Record `p` as an interior node of the segment `(q0, q1)` of string `gs`, if it
# lies strictly inside it. Exact, and endpoint-equal `p` is rejected first.
@inline function _record_vertex_on_segment!(m::Manifold, table::NodeTable, seg_nodes,
        p, q0, q1, gs::Int32, k::Int32; exact)
    (p == q0 || p == q1) && return nothing
    rk_point_on_segment(m, p, q0, q1; exact) || return nothing
    _record_interior!(seg_nodes, gs, k, _intern_node!(table, vertex_node(p)))
    return nothing
end

#=
Every unordered pair of distinct leaves of one spatial index whose extents
intersect, each visited once, `f(i1, i2)` with `i1 < i2` in leaf order.

`dual_depth_first_search(tree, tree)` would do the same job, but it visits every
unordered pair TWICE and every self-pair once, which the caller then has to
filter at the leaf — after the traversal has already paid for it. Splitting the
recursion into the within-child and cross-child halves removes that work instead
of discarding it, and is worth ~2x of the self-noding pass on dense real rings
(0.258 s -> 0.137 s over the 46 Natural Earth 10 m sides of the corpus quoted
above), where an input's own segments overlap each other constantly and there is
no pruning to be had.

`buf` materializes one leaf node's `(index, extent)` pairs so they can be walked
pairwise; it is threaded through the recursion and reused, so the whole traversal
allocates once. That makes `f` non-reentrant here — it must not start another
search — which neither caller does.
=#
function _self_pair_search(f::F, m::Manifold, tree, clip = nothing) where {F}
    buf = Tuple{Int, _segment_extent_type(m)}[]
    _self_pair_search(f, _clip_pair_pred(clip), tree, buf, clip)
    return nothing
end

# Whether an extent (a tree node's, or a leaf segment's) meets the clip box.
# Closed-interval, exactly like `_seg_in_clip` in split.jl — and on the same
# quantity, since a planar segment extent IS its endpoint bbox.
@inline _ext_in_clip(::Nothing, e) = true
@inline function _ext_in_clip(clip::Extents.Extent, e)
    (xlo, xhi) = clip.X
    (ylo, yhi) = clip.Y
    (ex_lo, ex_hi) = e.X
    (ey_lo, ey_hi) = e.Y
    return (ex_lo <= xhi) & (ex_hi >= xlo) & (ey_lo <= yhi) & (ey_hi >= ylo)
end

# The pair predicate: extents must meet each other AND both meet the clip box, so
# the same test prunes the cross-child `dual_depth_first_search` descent too.
_clip_pair_pred(::Nothing) = Extents.intersects
_clip_pair_pred(clip) = (ea, eb) ->
    Extents.intersects(ea, eb) && _ext_in_clip(clip, ea) && _ext_in_clip(clip, eb)

function _self_pair_search(f::F, pred::PR, node::N, buf::Vector, clip) where {F, PR, N}
    if SpatialTreeInterface.isleaf(node)
        empty!(buf)
        for ie in SpatialTreeInterface.child_indices_extents(node)
            #-- a pruned segment can pair with nothing that matters; drop it here
            #-- rather than in `pred`, so it does not cost a comparison per pair
            _ext_in_clip(clip, ie[2]) && push!(buf, ie)
        end
        @inbounds for a in 1:length(buf), b in (a + 1):length(buf)
            pred(buf[a][2], buf[b][2]) && f(buf[a][1], buf[b][1])
        end
    else
        nc = SpatialTreeInterface.nchild(node)
        for a in 1:nc
            ca = SpatialTreeInterface.getchild(node, a)
            ea = SpatialTreeInterface.node_extent(ca)
            #-- a subtree wholly outside the clip box holds only pruned segments
            _ext_in_clip(clip, ea) || continue
            #-- pairs inside this child ...
            _self_pair_search(f, pred, ca, buf, clip)
            #-- ... then pairs spanning it and a later sibling
            for b in (a + 1):nc
                cb = SpatialTreeInterface.getchild(node, b)
                pred(ea, SpatialTreeInterface.node_extent(cb)) &&
                    SpatialTreeInterface.dual_depth_first_search(f, pred, ca, cb)
            end
        end
    end
    return nothing
end

#=
## Rejecting a candidate pair in Float64 (design §2.3, performance)

`_classify_pair!` costs four adaptive-exact orientations, and on offset-curve
linework the overwhelming majority of candidates record nothing. Two shapes
dominate:

  * **segments sharing an endpoint** — every ring-adjacent pair, and every pair
    inside a fillet arc, whose boxes necessarily overlap;
  * **boxes that overlap while the segments stay apart** — a long inside-turn
    spoke whose thin diagonal box crosses the boxes of many short arc segments
    it never comes near.

`_pair_needs_kernel` rejects both in plain Float64. Every sign it acts on
carries Shewchuk's stage-A certificate — the same test `AdaptivePredicates`
applies before escalating, so a certified sign is the exact sign and an
uncertain one falls through to the kernel untouched. The recorded output is
therefore identical pair for pair, which `perf/noder` asserts on every
workload. The filter is planar and exact-mode only: `exact = False()` runs a
deliberately sloppier kernel that this must not second-guess.

Shared endpoints need the collinearity caveat. Two segments meeting at a shared
endpoint `q` record nothing *unless* they overlap collinearly — a hairpin, which
raw offset curves through input vertices really do contain. Otherwise the only
incidence flags set are the ones at `q`, and `_classify_pair!`'s own
`a0 != b0 && a0 != b1` guards already drop those; a second incidence would put
two points of one segment on the other, forcing collinearity. So one certified
orientation of `(other_a, q, other_b)` decides it: certainly non-zero means
certainly nothing to record. Four orientations become one.
=#

#-- Shewchuk's `ccwerrboundA` for Float64, (3 + 16ε)ε with ε = 2^-53.
const _ORIENT_ERRBOUND_A = (3.0 + 16.0 * (eps(Float64) / 2)) * (eps(Float64) / 2)

# `orient(a, b, c)` in plain Float64 with its stage-A certificate: `certain`
# means the returned value has the sign of the exact determinant (including a
# certified zero). This is verbatim the filter stage of
# `AdaptivePredicates.orient2`, so the two never disagree.
@inline function _orient_filter(ax, ay, bx, by, cx, cy)
    detleft = (ax - cx) * (by - cy)
    detright = (ay - cy) * (bx - cx)
    det = detleft - detright
    if detleft > 0
        detright <= 0 && return (det, true)
        detsum = detleft + detright
    elseif detleft < 0
        detright >= 0 && return (det, true)
        detsum = -detleft - detright
    else
        return (det, true)
    end
    errb = _ORIENT_ERRBOUND_A * detsum
    return (det, (det >= errb) | (-det >= errb))
end

# `p`, `r` are the two segments' free endpoints, `q` their shared one.
@inline function _shared_endpoint_needs_kernel(px, py, qx, qy, rx, ry)
    d, certain = _orient_filter(px, py, qx, qy, rx, ry)
    return !(certain & (d != 0))
end

@inline _pair_needs_kernel(::Manifold, ::BoolsAsTypes, a0, a1, b0, b1) = true

@inline function _pair_needs_kernel(::Planar, ::True, a0, a1, b0, b1)
    ax0, ay0 = GI.x(a0), GI.y(a0); ax1, ay1 = GI.x(a1), GI.y(a1)
    bx0, by0 = GI.x(b0), GI.y(b0); bx1, by1 = GI.x(b1), GI.y(b1)
    #-- (1) a shared endpoint: only a collinear overlap can record anything
    if (ax0 == bx0) & (ay0 == by0)
        return _shared_endpoint_needs_kernel(ax1, ay1, ax0, ay0, bx1, by1)
    elseif (ax0 == bx1) & (ay0 == by1)
        return _shared_endpoint_needs_kernel(ax1, ay1, ax0, ay0, bx0, by0)
    elseif (ax1 == bx0) & (ay1 == by0)
        return _shared_endpoint_needs_kernel(ax0, ay0, ax1, ay1, bx1, by1)
    elseif (ax1 == bx1) & (ay1 == by1)
        return _shared_endpoint_needs_kernel(ax0, ay0, ax1, ay1, bx0, by0)
    end
    #-- (2) one segment certainly strictly on one side of the other's line
    d0, c0 = _orient_filter(bx0, by0, bx1, by1, ax0, ay0)
    d1, c1 = _orient_filter(bx0, by0, bx1, by1, ax1, ay1)
    (c0 & c1 & (d0 != 0) & (d1 != 0) & ((d0 > 0) == (d1 > 0))) && return false
    e0, f0 = _orient_filter(ax0, ay0, ax1, ay1, bx0, by0)
    e1, f1 = _orient_filter(ax0, ay0, ax1, ay1, bx1, by1)
    (f0 & f1 & (e0 != 0) & (e1 != 0) & ((e0 > 0) == (e1 > 0))) && return false
    return true
end

# Function barrier: statically-typed classification of one candidate pair
# (the do-blocks above are dynamic closures over the tree traversal). `gsa`/`gsb`
# are the *global* string indices (into the arrangement's `segstrings`) of the
# two segments' parents; `sa`/`sb` index the lists actually passed in.
function _classify_pair!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ssa, sa::Int, gsa::Int32, ksa::Int32,
        ssb, sb::Int, gsb::Int32, ksb::Int32; exact) where {P}
    a0 = ssa[sa].pts[ksa]; a1 = ssa[sa].pts[ksa + 1]
    b0 = ssb[sb].pts[ksb]; b1 = ssb[sb].pts[ksb + 1]

    _pair_needs_kernel(m, booltype(exact), a0, a1, b0, b1) || return nothing
    cls = rk_classify_intersection(m, a0, a1, b0, b1; exact)
    kind = cls.kind
    if kind == SS_DISJOINT
        return nothing
    elseif kind == SS_PROPER
        #-- a proper crossing is strictly interior to BOTH segments
        nid = _intern_node!(table, crossing_node(a0, a1, b0, b1))
        _record_interior!(seg_nodes, gsa, ksa, nid)
        _record_interior!(seg_nodes, gsb, ksb, nid)
        return nothing
    end
    #-- SS_TOUCH / SS_COLLINEAR: every intersection point is an input vertex,
    #-- reported via the incidence flags (design §2.3, S1 census). The claim is
    #-- load-bearing; assert it.
    @assert (cls.a0_on_b || cls.a1_on_b || cls.b0_on_a || cls.b1_on_a) """
        touch/collinear intersection with no vertex incidence flag — \
        the touch=input-vertex invariant (design §2.3) is violated"""
    #-- a?/b? are vertices of A/B; record on the OTHER segment only where the
    #-- vertex lies strictly in that segment's interior (not at its endpoints).
    if cls.a0_on_b && a0 != b0 && a0 != b1
        _record_interior!(seg_nodes, gsb, ksb, _intern_node!(table, vertex_node(a0)))
    end
    if cls.a1_on_b && a1 != b0 && a1 != b1
        _record_interior!(seg_nodes, gsb, ksb, _intern_node!(table, vertex_node(a1)))
    end
    if cls.b0_on_a && b0 != a0 && b0 != a1
        _record_interior!(seg_nodes, gsa, ksa, _intern_node!(table, vertex_node(b0)))
    end
    if cls.b1_on_a && b1 != a0 && b1 != a1
        _record_interior!(seg_nodes, gsa, ksa, _intern_node!(table, vertex_node(b1)))
    end
    return nothing
end
