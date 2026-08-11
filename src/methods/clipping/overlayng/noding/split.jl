# # Stages 2 + 4 — order and split (design §2.5, §2.1)
#
# For every parent segment of every string, order its interior nodes along the
# segment (stage 2, `rk_compare_along_segment`), cap the chain with the segment's
# endpoint vertex nodes, and emit one `NodedEdge` per non-degenerate link
# (stage 4). Interior node lists are already deduped by `node_identity.jl`; the
# consecutive-dedup here guards the endpoint↔interior joins.
#
# Ordering is zero-cost for the <2-node case (the overwhelming majority — S1: no
# real segment acquired ≥2 interior nodes), which is the only case that touches
# the kernel ordering predicate.

#=
## Clip pruning (the construct-free `RingClipper`)

SCOPE, because it is easy to over-read this section. Clip pruning is inert
outside a narrow case, and any behaviour observed there is NOT its doing:

  * it is PLANAR-ONLY — `_overlay_envelopes` returns `(nothing, nothing)` for
    every non-`Planar` manifold, so no spherical query ever receives a box; and
  * it applies to INTERSECTION and DIFFERENCE only — `_overlay_clip_envelopes`
    returns `(nothing, nothing)` for UNION and SYMDIFFERENCE, which need every
    edge of both inputs.

So both `clip_a` and `clip_b` are `nothing` for a spherical op, and for a union
on either manifold, and every `_seg_in_clip` / `_ext_in_clip` / `_pt_in_clip`
test below degenerates to its `::Nothing` method. In particular, nothing in this
file participates in a spherical cascaded union.

`clip_a` / `clip_b` are optional per-side bounding boxes. A parent segment whose
own bounding box misses its side's box emits NO `NodedEdge` at all — it is not
clipped, not shortened, and no coordinate is synthesized on the box, so nothing
constructed can reach a decision (design §0). The kept sub-chains of a ring are
simply left OPEN: a ring that leaves and re-enters the box becomes one or more
disconnected chains, whose ends are ordinary degree-1 graph nodes.

Two facts make that sound, and both are properties of *this* engine rather than
of JTS/GEOS (where the same trick would be wrong):

  * the winding authority is per segment STRING, not per edge —
    `EdgeSourceInfo.depth_delta` is computed once at ingest from the original
    ring (`_ring_material_interior_on_left`), so removing edges cannot perturb
    it; and
  * an edge whose location is left unknown falls back to point-in-area against
    the ORIGINAL input (`_label_disconnected_edges!` through `_OverlayInput`,
    design §3 amendment 7), never against the surviving linework.

The caller's per-op envelope rule (`_overlay_clip_envelopes`) supplies the third
fact — that every in-result edge lies inside the box — and with it the invariant
this stage exports:

    every node lying inside the clip box keeps its FULL star.

A pruned segment's bbox misses the box, so every point of it (its endpoints
included) is outside; contrapositively, every segment incident on a point in the
box has a bbox meeting the box and survives. Only a node OUTSIDE the box can
therefore have lost an edge, and `truncated` marks exactly those, so the one pass
that reads a star as a closed cycle (`_propagate_area_locations!`) can skip them.

`_truncated_bits` reads that off the node keys directly: a node is outside the
box iff its key is a VERTEX key whose coordinate is outside. Crossing keys are
never marked, and need not be — a crossing lies on a segment of each input, hence
inside both envelopes, hence inside the box. No coordinate is constructed to
decide this: a vertex key stores the original input vertex verbatim, so the test
reads an input coordinate, not an emitted one. (An earlier version tracked the
pruned segments' own endpoint ids instead. That is a strictly smaller set, but it
has to be maintained by every stage that can drop a node — which the self-noding
prune in collect.jl then broke — and it cost a hash per pruned endpoint. The
positional rule is a superset, is O(#nodes) of pointer-free comparisons, and
cannot go stale.)
=#

@inline _seg_in_clip(::Nothing, p0, p1) = true

@inline function _seg_in_clip(clip::Extents.Extent, p0::NTuple{2, Float64},
        p1::NTuple{2, Float64})
    (xlo, xhi) = clip.X
    (ylo, yhi) = clip.Y
    #-- CLOSED-interval overlap: a segment merely touching the box is kept
    return (min(p0[1], p1[1]) <= xhi) & (max(p0[1], p1[1]) >= xlo) &
           (min(p0[2], p1[2]) <= yhi) & (max(p0[2], p1[2]) >= ylo)
end

function _split_edges!(m::Manifold, table::NodeTable{P}, seg_nodes,
        segstrings::Vector{RelateSegmentString{P}}, na::Int32,
        clip_a, clip_b; exact) where {P}
    edges = NodedEdge[]
    chain = Int32[]
    for gsi in eachindex(segstrings)
        ss = segstrings[gsi]
        pts = ss.pts
        gsi32 = Int32(gsi)
        clip = gsi32 <= na ? clip_a : clip_b
        for k in 1:(length(pts) - 1)
            #-- clip prune: no edge, and no endpoint interning, for this segment
            _seg_in_clip(clip, pts[k], pts[k + 1]) || continue
            interior = get(seg_nodes, (gsi32, Int32(k)), nothing)
            if interior !== nothing && length(interior) >= 2
                _order_along_segment!(m, interior, pts[k], pts[k + 1], table; exact)
            end
            #-- build the node chain: start vertex, ordered interior, end vertex
            lo = _intern_node!(table, vertex_node(pts[k]))
            hi = _intern_node!(table, vertex_node(pts[k + 1]))
            empty!(chain)
            push!(chain, lo)
            if interior !== nothing
                append!(chain, interior)
            end
            push!(chain, hi)
            #-- emit one edge per link, skipping any zero-length (coincident) join
            for c in 1:(length(chain) - 1)
                a = chain[c]; b = chain[c + 1]
                a == b && continue
                push!(edges, NodedEdge(gsi32, Int32(k), a, b))
            end
        end
    end
    return edges, _truncated_bits(table, clip_a, clip_b)
end

@inline _pt_in_clip(::Nothing, p) = true
@inline function _pt_in_clip(clip::Extents.Extent, p::NTuple{2, Float64})
    (xlo, xhi) = clip.X
    (ylo, yhi) = clip.Y
    return (p[1] >= xlo) & (p[1] <= xhi) & (p[2] >= ylo) & (p[2] <= yhi)
end

# The nodes whose star clip pruning can have thinned: those outside a clip box.
# Empty (and untouched) whenever no box was given, so the unclipped arrangement
# carries no flags and the labeller's guard stays inert.
_truncated_bits(::NodeTable, ::Nothing, ::Nothing) = BitVector()
function _truncated_bits(table::NodeTable, clip_a, clip_b)
    keys = table.keys
    truncated = falses(length(keys))
    any_set = false
    @inbounds for i in eachindex(keys)
        k = keys[i]
        #-- a crossing lies on a segment of each input, so it is inside both
        #-- envelopes and inside the box; only vertex keys can be outside
        k.is_crossing && continue
        if !(_pt_in_clip(clip_a, k.pt) && _pt_in_clip(clip_b, k.pt))
            truncated[i] = true
            any_set = true
        end
    end
    return any_set ? truncated : BitVector()
end

# Sort a segment's interior node ids along the oriented segment (s0, s1). The
# kernel predicate is the sole authority; by construction (stage 3 ran first) no
# two distinct ids coincide, so a `0` comparison is an invariant violation.
function _order_along_segment!(m::Manifold, ids::Vector{Int32}, s0, s1,
        table::NodeTable; exact)
    sort!(ids; lt = (a, b) ->
        rk_compare_along_segment(m, s0, s1, table.keys[a], table.keys[b]; exact) < 0)
    #-- assert no residual coincidence among ordered neighbours (design §2.5)
    @assert begin
        ok = true
        for c in 2:length(ids)
            if rk_compare_along_segment(m, s0, s1, table.keys[ids[c - 1]], table.keys[ids[c]]; exact) == 0
                ok = false; break
            end
        end
        ok
    end "two distinct node ids coincide along a segment after node-identity merging (design §2.5)"
    return nothing
end
