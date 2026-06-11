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
[`rk_classify_intersection`](@ref) and add a [`NodeSection`](@ref) pair (one
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
