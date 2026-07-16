# # Half-edge substrate — the winged-edge algebra (port of JTS `edgegraph/HalfEdge`)
#
# Phase 2a of the OverlayNG port (design doc §3). A port of
# `edgegraph/HalfEdge.java`, adapted to integer-indexed storage: half-edges live
# in a `Vector` and refer to each other by 1-based index (`0` = null), not by
# object reference — the spike-validated representation ("indices, not object
# references, proved painless"). Directed half-edges occur in symmetric pairs;
# each stores the index of its symmetric partner (`sym`) and of the next
# half-edge CCW around its origin's star.
#
# Field contract (every element of the edge store must provide these mutable
# fields; `OverlayEdge` in overlay_graph.jl is the sole concrete implementation):
#   - `origin :: Int32`  — the origin node id (into the arrangement's node table).
#   - `sym    :: Int32`  — index of the symmetric (oppositely-oriented) half-edge.
#   - `o_next :: Int32`  — index of the next half-edge CCW around the origin, with
#                          the same origin (JTS `oNext`); a degree-1 origin points
#                          it at the half-edge itself.
#   - `dir_pt :: P`      — the direction point for angular ordering: the parent
#                          segment's far endpoint, an original kernel vertex
#                          (design §3 amendment 1) — never a constructed coordinate.
#
# Adaptation note (winged-edge field choice): JTS stores `next` (the next edge CCW
# around the DESTINATION) and derives `oNext() == sym.next`. We store `o_next`
# (around the ORIGIN) directly, because it is the primitive every overlay
# traversal uses and it makes `oNext`/`degree`/`prev` index-only with no `sym`
# indirection; JTS's `next` is recoverable as `next(e) == o_next(sym(e))`. We
# bulk-build each node's star and sort it once with the exact angular comparator
# (the spike's `_order_stars!`), in place of JTS's incremental `insert`: the exact
# comparator is costly, so one `O(d log d)` sort per node beats `insert`'s repeated
# rescans, and integer storage makes a bulk sort trivial.
#
# Everything here is internal to GeometryOps — nothing is exported.

# ## Core edge algebra (integer-indexed ports of the HalfEdge methods)

# The symmetric (oppositely-oriented) partner (JTS `sym`).
@inline he_sym(edges, i::Integer) = @inbounds edges[i].sym

# The next half-edge CCW around the origin, same origin (JTS `oNext`).
@inline he_onext(edges, i::Integer) = @inbounds edges[i].o_next

# The origin node id.
@inline he_origin(edges, i::Integer) = @inbounds edges[i].origin

# The destination node id: the origin of the sym edge (JTS `dest`).
@inline he_dest(edges, i::Integer) = @inbounds edges[he_sym(edges, i)].origin

# The next half-edge CCW around the destination (JTS `next`), recovered from the
# stored `o_next`: `next(e) == oNext(sym(e))`.
@inline he_next(edges, i::Integer) = he_onext(edges, he_sym(edges, i))

# The previous half-edge CW around the origin, with that vertex as its
# destination. Always `he_next(he_prev(e)) == e`. Port of JTS `prev` (scan the
# origin star, return the last edge's sym).
function he_prev(edges, i::Integer)
    curr = Int32(i); prev = Int32(i)
    while true
        prev = curr
        curr = he_onext(edges, curr)
        curr == i && break
    end
    return he_sym(edges, prev)
end

# The degree of the origin vertex: the number of half-edges originating there
# (port of JTS `degree`).
function he_degree(edges, i::Integer)
    d = 0
    e = Int32(i)
    while true
        d += 1
        e = he_onext(edges, e)
        e == i && break
    end
    return d
end

# The first node (degree != 2) reached walking `prev` from this edge, or `0` if
# the edge is part of a ring with no such node (port of JTS `prevNode`).
function he_prev_node(edges, i::Integer)
    e = Int32(i)
    while he_degree(edges, e) == 2
        e = he_prev(edges, e)
        e == i && return Int32(0)
    end
    return e
end

# The half-edge originating at this edge's origin with the given destination node
# id, or `0` if none (port of JTS `find`, on node ids rather than coordinates).
function he_find(edges, i::Integer, dest::Integer)
    o = Int32(i)
    while true
        he_dest(edges, o) == dest && return o
        o = he_onext(edges, o)
        o == i && break
    end
    return Int32(0)
end

# ## Symmetric-pair linkage and star ordering

# Link a freshly created symmetric pair `(i, j)`: set the sym pointers and the
# single-segment star (each edge's `o_next` is itself, so a degree-1 origin's
# `oNext` is the edge itself — JTS `link`, restated for the origin star). Both are
# overwritten for interior-star edges once their origin's star is ordered.
@inline function he_link!(edges, i::Integer, j::Integer)
    @inbounds edges[i].sym = Int32(j)
    @inbounds edges[j].sym = Int32(i)
    @inbounds edges[i].o_next = Int32(i)
    @inbounds edges[j].o_next = Int32(j)
    return nothing
end

# Angular comparison of two half-edges that share an origin, about that origin's
# symbolic apex (the node's `NodeKey`), via the exact kernel comparator
# (`rk_compare_edge_dir`, design §3 amendment 1). Foreign directions at a
# coincidence-merged apex take the kernel's exact-rational slow path — no special
# casing here.
@inline function he_compare_angular(m::Manifold, edges, keys, i::Integer, j::Integer; exact)
    apex = @inbounds keys[he_origin(edges, i)]
    return rk_compare_edge_dir(m, apex, (@inbounds edges[i].dir_pt), (@inbounds edges[j].dir_pt); exact)
end

# Order one node's star of outgoing half-edge indices CCW and wire the origin ring
# so that `he_onext` walks it (the spike's `_order_stars!`). `star` is mutated
# into CCW order. All members must share the origin whose apex is `keys[origin]`.
function he_order_star!(m::Manifold, edges, keys, star::Vector{Int32}; exact)
    n = length(star)
    if n == 1
        @inbounds edges[star[1]].o_next = star[1]
        return nothing
    end
    sort!(star; lt = (i, j) -> he_compare_angular(m, edges, keys, i, j; exact) < 0)
    @inbounds for t in 1:n
        edges[star[t]].o_next = star[t == n ? 1 : t + 1]
    end
    return nothing
end
