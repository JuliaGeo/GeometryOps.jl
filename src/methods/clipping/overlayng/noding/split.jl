# # Stages 2 + 4 — order and split (design §2.5, §2.1)
#
# For every parent segment of every string, order its interior nodes along the
# segment (stage 2, `rk_compare_along_segment`), cap the chain with the segment's
# endpoint vertex nodes, and emit one [`NodedEdge`](@ref) per non-degenerate link
# (stage 4). Interior node lists are already deduped by `node_identity.jl`; the
# consecutive-dedup here guards the endpoint↔interior joins.
#
# Ordering is zero-cost for the <2-node case (the overwhelming majority — S1: no
# real segment acquired ≥2 interior nodes), which is the only case that touches
# the kernel ordering predicate.

function _split_edges!(m::Manifold, table::NodeTable{P}, seg_nodes,
        segstrings::Vector{RelateSegmentString{P}}; exact) where {P}
    edges = NodedEdge[]
    chain = Int32[]
    for gsi in eachindex(segstrings)
        ss = segstrings[gsi]
        pts = ss.pts
        gsi32 = Int32(gsi)
        for k in 1:(length(pts) - 1)
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
    return edges
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
