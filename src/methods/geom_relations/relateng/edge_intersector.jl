# # RelateNG edge segment intersector
#
# Port of JTS `EdgeSegmentIntersector.java`: tests segments of
# [`RelateSegmentString`](@ref)s and, if they intersect, adds the
# intersection(s) to the [`TopologyComputer`](@ref).
#
# The Java class wraps a `RobustLineIntersector`, which *constructs* the
# intersection coordinate(s) and then filters endpoint intersections through
# `RelateSegmentString.isContainingSegment` so each one is processed once
# only, on its canonical segments. Here the construction is replaced by the
# symbolic classification [`rk_classify_intersection`](@ref) (design D2):
#
# - `SS_PROPER` ↔ `li.isProper()`: the node is identified by its
#   [`crossing_node`](@ref) key (the canonicalized defining segment pair) —
#   no coordinate is ever computed. Proper intersections lie on a unique
#   segment pair, so (as in Java) no canonicality check is needed.
# - `SS_TOUCH`/`SS_COLLINEAR` ↔ the non-proper point and collinear outcomes:
#   every JTS intersection point is an *input vertex* lying on the other
#   segment, which is exactly what the `SegSegClass` incidence flags
#   (`a0_on_b`, `a1_on_b`, `b0_on_a`, `b1_on_a`) report. Each true flag is a
#   candidate intersection point; flags naming the same coordinate (e.g. a
#   shared endpoint `a0 == b0` sets both `a0_on_b` and `b0_on_a`) denote ONE
#   geometric point. A collinear overlap has up to two distinct points (the
#   ends of the overlap interval), matching `li.getIntersectionNum() == 2`.
#
# Per distinct intersection point the Java loop adds one `NodeSection` on
# ssA AND one on ssB (`topologyComputer.addIntersection(nsa, nsb)`); the
# point is processed only if BOTH strings contain it canonically.
#
# The Java `isDone()` (the noder's early-exit hook) is `is_result_known` on
# the computer; the Task 20 enumerator consults it directly.
#
# Method order parallels the Java file, so this file diffs against its
# counterpart.

"""
    process_intersections!(tc::TopologyComputer, ss0, seg_index0, ss1, seg_index1; m, exact)

Classify the intersection of segment `seg_index0` of `ss0` with segment
`seg_index1` of `ss1` and record any intersections in `tc`. The strings are
ordered so the A geometry's string is processed first (the computer's A/B
matrix updates rely on it).

Port of EdgeSegmentIntersector.processIntersections.
"""
function process_intersections!(tc::TopologyComputer,
        ss0::RelateSegmentString, seg_index0::Integer,
        ss1::RelateSegmentString, seg_index1::Integer;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    #-- don't intersect a segment with itself
    ss0 === ss1 && seg_index0 == seg_index1 && return nothing

    if ss0.is_a
        add_intersections!(tc, ss0, seg_index0, ss1, seg_index1; m, exact)
    else
        add_intersections!(tc, ss1, seg_index1, ss0, seg_index0; m, exact)
    end
    return nothing
end

"""
    add_intersections!(tc::TopologyComputer, ssA, seg_index_a, ssB, seg_index_b; m, exact)

Classify the intersection of one segment pair via
`rk_classify_intersection` and add a [`NodeSection`](@ref) pair (one
on `ssA`, one on `ssB`) to `tc` for each distinct intersection point:

- `SS_DISJOINT`: nothing to add.
- `SS_PROPER`: one section pair at the symbolic [`crossing_node`](@ref).
- `SS_TOUCH`/`SS_COLLINEAR`: a section pair at the [`vertex_node`](@ref) of
  each distinct flagged vertex — but only when both segments contain the
  vertex canonically ([`_is_canonical_incidence`](@ref)), which ensures
  endpoint intersections are added once only across adjacent segments.

Port of EdgeSegmentIntersector.addIntersections (private).
"""
function add_intersections!(tc::TopologyComputer,
        ssA::RelateSegmentString, seg_index_a::Integer,
        ssB::RelateSegmentString, seg_index_b::Integer;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    a0 = ssA.pts[seg_index_a]
    a1 = ssA.pts[seg_index_a + 1]
    b0 = ssB.pts[seg_index_b]
    b1 = ssB.pts[seg_index_b + 1]

    cls = rk_classify_intersection(m, a0, a1, b0, b1; exact)

    cls.kind == SS_DISJOINT && return nothing

    if cls.kind == SS_PROPER
        #-- a proper intersection lies on a unique segment pair, so it needs
        #-- no canonicality check (and its node is purely symbolic)
        node = crossing_node(a0, a1, b0, b1)
        nsa = create_node_section(ssA, seg_index_a, node)
        nsb = create_node_section(ssB, seg_index_b, node)
        add_intersection!(tc, nsa, nsb)
        return nothing
    end

    #-- SS_TOUCH / SS_COLLINEAR: each flagged vertex is a candidate
    #-- intersection point. Flags naming the same coordinate (a shared
    #-- endpoint flags both strings) are deduped via the normalized vertex
    #-- key; at most two distinct points exist (a collinear overlap's ends).
    seen1 = seen2 = nothing
    for (flag, pt) in ((cls.a0_on_b, a0), (cls.a1_on_b, a1),
                       (cls.b0_on_a, b0), (cls.b1_on_a, b1))
        flag || continue
        node = vertex_node(pt)
        (node === seen1 || node === seen2) && continue
        if seen1 === nothing
            seen1 = node
        else
            seen2 = node
        end
        #-- ensure endpoint intersections are added once only, for their
        #-- canonical segments
        if _is_canonical_incidence(ssA, seg_index_a, pt) &&
                _is_canonical_incidence(ssB, seg_index_b, pt)
            nsa = create_node_section(ssA, seg_index_a, node)
            nsb = create_node_section(ssB, seg_index_b, node)
            add_intersection!(tc, nsa, nsb)
        end
    end
    return nothing
end

"""
    _is_canonical_incidence(ss::RelateSegmentString, seg_index::Integer, pt)

The once-only rule for vertex incidences: whether segment `seg_index` of
`ss` is the canonical owner of the intersection point `pt`. Segments are
half-closed — a segment owns its start vertex but not its end vertex, except
for the final segment of a non-closed string, which also owns its endpoint;
in a closed ring the closing vertex is owned by the first segment
(wraparound). This attributes every vertex of a segment string to exactly
one of its segments, so an endpoint intersection enumerated against both
incident segments produces node sections once, not twice.

Encodes the canonicality semantics of the Java
`addIntersections`/`RelateSegmentString.isContainingSegment` pairing;
delegates to [`is_containing_segment`](@ref) (its direct port).
"""
_is_canonical_incidence(ss::RelateSegmentString, seg_index::Integer, pt) =
    is_containing_segment(ss, seg_index, pt)

#==========================================================================
## Edge set enumeration

Replaces JTS `EdgeSetIntersector.java` (HPRtree + monotone chains) and the
mutual A x B pruning of `MCIndexSegmentSetMutualIntersector`: enumerate
every segment pair (one segment from an A string, one from a B string)
whose extents interact and feed it through `process_intersections!`.

The accelerator strategy mirrors the clipping pattern
(`foreach_pair_of_maybe_intersecting_edges_in_order` in
clipping_processor.jl), reusing its `IntersectionAccelerator` types and the
`GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS` size threshold.

The Java noder's `isDone()` early-exit hook lands here: after each
processed pair the enumerator consults `is_result_known(computer)` and
stops as soon as the predicate value is determined. On the tree path the
traversal is terminated by returning `Action(:full_return, nothing)` from
the callback — `dual_depth_first_search` processes `LoopStateMachine`
actions via `@controlflow`, and `:full_return` propagates out of the whole
recursion (a plain `:break` would only exit the innermost leaf loop), so
no exception is needed.
==========================================================================#

"""
    process_edge_intersections!(computer, ssa_list, ssb_list,
        accelerator = AutoAccelerator(); m, exact)

Enumerate all extent-interacting segment pairs between the A-side segment
strings `ssa_list` and the B-side segment strings `ssb_list`, feeding each
pair through [`process_intersections!`](@ref) so the intersections are
recorded on `computer`.

`accelerator` selects the enumeration strategy:

- `NestedLoop`: a plain double loop over
  string pairs and segment pairs, with a per-pair segment-extent
  disjointness skip (on `Planar`).
- Any tree-backed accelerator (e.g. `DoubleSTRtree`): a spatial index
  (`_relate_edge_index`, currently a `NaturalIndex`) is built over the
  per-segment extents of each side and traversed with
  `SpatialTreeInterface.dual_depth_first_search`
  under the `Extents.intersects` predicate.
- [`AutoAccelerator`](@ref): picks `NestedLoop` below the clipping size
  threshold (`GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS`) and on
  non-`Planar` manifolds (planar extent trees are not valid there), and the
  tree path otherwise.

After each processed pair `is_result_known(computer)` is consulted and the
enumeration stops early once the predicate value is determined (the port of
the Java noder's `isDone()` hook used by `EdgeSetIntersector.process`).

!!! warning
    This enumerates A×B pairs only. JTS's `EdgeSetIntersector` also feeds
    A×A and B×B pairs (self-noding) with an id-ordering guard so each
    unordered pair is processed once. Calling this with the same list on
    both sides would process every pair twice — the engine's
    `computeAtEdges` port uses [`process_self_intersections!`](@ref) for
    the self-pair path instead.
"""
process_edge_intersections!(tc::TopologyComputer,
        ssa_list::AbstractVector{<:RelateSegmentString},
        ssb_list::AbstractVector{<:RelateSegmentString};
        m::Manifold = _manifold(tc), exact = _exact(tc)) =
    process_edge_intersections!(tc, ssa_list, ssb_list, AutoAccelerator(); m, exact)

# AutoAccelerator: pick the strategy from the manifold and the segment
# counts, mirroring the clipping selection
# (`foreach_pair_of_maybe_intersecting_edges_in_order(::AutoAccelerator)`).
function process_edge_intersections!(tc::TopologyComputer,
        ssa_list::AbstractVector{<:RelateSegmentString},
        ssb_list::AbstractVector{<:RelateSegmentString},
        ::AutoAccelerator;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    return process_edge_intersections!(tc, ssa_list, ssb_list,
        _select_edge_set_accelerator(m, ssa_list, ssb_list); m, exact)
end

# STRtrees over planar extents are only valid on the Planar manifold, and
# below the clipping threshold the nested loop wins anyway.
function _select_edge_set_accelerator(::Planar, ssa_list, ssb_list)
    na = _total_segment_count(ssa_list)
    nb = _total_segment_count(ssb_list)
    if na < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS &&
            nb < GEOMETRYOPS_NO_OPTIMIZE_EDGEINTERSECT_NUMVERTS
        return NestedLoop()
    else
        return DoubleSTRtree()
    end
end
_select_edge_set_accelerator(::Manifold, ssa_list, ssb_list) = NestedLoop()

_total_segment_count(ss_list) =
    sum(ss -> length(ss.pts) - 1, ss_list; init = 0)

# NestedLoop path: double loop over string pairs x segment pairs, skipping
# pairs whose segment extents are disjoint (the pruning that the monotone
# chains of Java's MCIndexSegmentSetMutualIntersector provide).
function process_edge_intersections!(tc::TopologyComputer,
        ssa_list::AbstractVector{<:RelateSegmentString},
        ssb_list::AbstractVector{<:RelateSegmentString},
        ::NestedLoop;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    for ssa in ssa_list, ssb in ssb_list
        for ia in 1:(length(ssa.pts) - 1)
            a0 = ssa.pts[ia]
            a1 = ssa.pts[ia + 1]
            for ib in 1:(length(ssb.pts) - 1)
                b0 = ssb.pts[ib]
                b1 = ssb.pts[ib + 1]
                _segment_envs_disjoint(m, a0, a1, b0, b1) && continue
                process_intersections!(tc, ssa, ia, ssb, ib; m, exact)
                #-- the Java noder's isDone() early-exit hook
                is_result_known(tc) && return nothing
            end
        end
    end
    return nothing
end

# Per-pair extent pruning is a planar coordinate comparison; on other
# manifolds it could wrongly discard interacting pairs (e.g. across the
# antimeridian), so prune nothing there.
_segment_envs_disjoint(::Planar, a0, a1, b0, b1) =
    min(a0[1], a1[1]) > max(b0[1], b1[1]) ||
    max(a0[1], a1[1]) < min(b0[1], b1[1]) ||
    min(a0[2], a1[2]) > max(b0[2], b1[2]) ||
    max(a0[2], a1[2]) < min(b0[2], b1[2])
_segment_envs_disjoint(::Manifold, a0, a1, b0, b1) = false

# The spatial index built over per-segment extents for the tree-accelerated
# paths (here and in the prepared mode of relate_ng.jl). A `NaturalIndex`
# rather than an `STRtree`: segments arrive in ring/line order, which is
# already spatially coherent, so the no-sort natural index (pure in-order
# hierarchical extent reduction) builds much faster while pruning the dual
# traversal almost as well. Both implement SpatialTreeInterface, so this is
# the only line to change to swap index structures.
_relate_edge_index(extents::Vector{<:Extents.Extent}) =
    NaturalIndex(extents; nodecapacity = 16)

# Tree path (any other accelerator, canonically DoubleSTRtree): a spatial
# index over the per-segment extents of each side, traversed simultaneously.
function process_edge_intersections!(tc::TopologyComputer,
        ssa_list::AbstractVector{<:RelateSegmentString},
        ssb_list::AbstractVector{<:RelateSegmentString},
        ::IntersectionAccelerator;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    extents_a, owners_a = _segment_extent_table(ssa_list)
    extents_b, owners_b = _segment_extent_table(ssb_list)
    (isempty(extents_a) || isempty(extents_b)) && return nothing
    tree_a = _relate_edge_index(extents_a)
    tree_b = _relate_edge_index(extents_b)
    SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree_a, tree_b) do ia, ib
        (sa, ka) = owners_a[ia]
        (sb, kb) = owners_b[ib]
        process_intersections!(tc, ssa_list[sa], ka, ssb_list[sb], kb; m, exact)
        #-- the Java noder's isDone() early-exit hook; :full_return
        #-- propagates out of the whole dual traversal via @controlflow
        is_result_known(tc) && return Action(:full_return, nothing)
        return nothing
    end
    return nothing
end

#==========================================================================
## Self-pair enumeration (the A×A / B×B side of JTS EdgeSetIntersector)

When `is_self_noding_required(tc)` holds, JTS's `computeEdgesAll` puts the
edges of BOTH inputs into one `EdgeSetIntersector`, whose `process` visits
every unordered pair of distinct monotone chains exactly once (the
`testChain.getId() <= queryChain.getId()` guard) — i.e. all A×B pairs plus
the A×A and B×B self pairs. The mutual A×B pairs are handled by
`process_edge_intersections!` above; `process_self_intersections!` is the
guarded self-pair path for one side's list.

The id-ordering guard becomes: unordered string pairs `si <= sj`, and
within a single string (`si == sj`) unordered segment pairs `ka < kb`
(never a segment with itself). Note JTS never compares a chain with
itself — safe there because a monotone chain cannot self-intersect; a
whole segment string can, so same-string segment pairs ARE enumerated
here. Trivial adjacent-segment endpoint touches are filtered by the same
canonical-incidence rule as in Java (`is_containing_segment`), so this
produces exactly the node sections the Java chain enumeration does.
==========================================================================#

"""
    process_self_intersections!(computer, ss_list,
        accelerator = AutoAccelerator(); m, exact)

Enumerate all extent-interacting segment pairs *within* the segment-string
list `ss_list` (each unordered pair once, never a segment against itself),
feeding each through [`process_intersections!`](@ref) — the self-noding
(A×A or B×B) counterpart of [`process_edge_intersections!`](@ref).

After each processed pair `is_result_known(computer)` is consulted and the
enumeration stops early once the predicate value is determined.
"""
process_self_intersections!(tc::TopologyComputer,
        ss_list::AbstractVector{<:RelateSegmentString};
        m::Manifold = _manifold(tc), exact = _exact(tc)) =
    process_self_intersections!(tc, ss_list, AutoAccelerator(); m, exact)

function process_self_intersections!(tc::TopologyComputer,
        ss_list::AbstractVector{<:RelateSegmentString},
        ::AutoAccelerator;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    return process_self_intersections!(tc, ss_list,
        _select_edge_set_accelerator(m, ss_list, ss_list); m, exact)
end

# NestedLoop path: unordered string pairs si <= sj; within one string the
# segment pairs are also unordered (ka < kb).
function process_self_intersections!(tc::TopologyComputer,
        ss_list::AbstractVector{<:RelateSegmentString},
        ::NestedLoop;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    for si in eachindex(ss_list)
        ssa = ss_list[si]
        for sj in si:lastindex(ss_list)
            ssb = ss_list[sj]
            for ia in 1:(length(ssa.pts) - 1)
                a0 = ssa.pts[ia]
                a1 = ssa.pts[ia + 1]
                ib0 = si == sj ? ia + 1 : 1
                for ib in ib0:(length(ssb.pts) - 1)
                    b0 = ssb.pts[ib]
                    b1 = ssb.pts[ib + 1]
                    _segment_envs_disjoint(m, a0, a1, b0, b1) && continue
                    process_intersections!(tc, ssa, ia, ssb, ib; m, exact)
                    #-- the Java noder's isDone() early-exit hook
                    is_result_known(tc) && return nothing
                end
            end
        end
    end
    return nothing
end

# Tree path: one index over the per-segment extents, dual-traversed with
# itself; the flat-index ordering `ia < ib` is the unordered-pair guard
# (excluding a segment against itself).
function process_self_intersections!(tc::TopologyComputer,
        ss_list::AbstractVector{<:RelateSegmentString},
        ::IntersectionAccelerator;
        m::Manifold = _manifold(tc), exact = _exact(tc))
    extents, owners = _segment_extent_table(ss_list)
    isempty(extents) && return nothing
    tree = _relate_edge_index(extents)
    SpatialTreeInterface.dual_depth_first_search(Extents.intersects, tree, tree) do ia, ib
        ia < ib || return nothing
        (sa, ka) = owners[ia]
        (sb, kb) = owners[ib]
        process_intersections!(tc, ss_list[sa], ka, ss_list[sb], kb; m, exact)
        #-- the Java noder's isDone() early-exit hook; :full_return
        #-- propagates out of the whole dual traversal via @controlflow
        is_result_known(tc) && return Action(:full_return, nothing)
        return nothing
    end
    return nothing
end

# Flat per-segment extent list for a segment-string list, with the offset
# table mapping each flat index back to (string index, segment index).
function _segment_extent_table(ss_list)
    extents = Extents.Extent{(:X, :Y), NTuple{2, NTuple{2, Float64}}}[]
    owners = NTuple{2, Int}[]
    nseg = _total_segment_count(ss_list)
    sizehint!(extents, nseg)
    sizehint!(owners, nseg)
    for (si, ss) in enumerate(ss_list)
        _push_segment_extents!(extents, owners, si, ss.pts)
    end
    return extents, owners
end

# Function barrier: dispatch once per segment string, so the per-segment
# loop stays statically typed even if `ss_list` has a non-concrete eltype.
function _push_segment_extents!(extents::Vector, owners::Vector, si::Int, pts::Vector)
    for k in 1:(length(pts) - 1)
        p = pts[k]
        q = pts[k + 1]
        push!(extents, Extents.Extent(X = minmax(p[1], q[1]), Y = minmax(p[2], q[2])))
        push!(owners, (si, k))
    end
    return nothing
end
