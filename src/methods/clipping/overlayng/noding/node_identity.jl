# # Stage 3 — node identity (design §2.4)
#
# Two tiers, no canonical key. Tier 1 (egal `NodeKey` equality) already ran
# during collect via `_intern_node!`. Tier 2 finds *geometric* coincidence across
# distinct keys — two segment pairs crossing at one point, or a crossing landing
# on a third string's vertex — which egal equality misses. A throwaway Float64
# proximity sweep flags candidate pairs; the exact `rk_nodes_coincide` confirms;
# a union-find over confirmed pairs (bounded by candidate count — zero on all
# real data measured, a handful on constructed degree-≥3 nodes) merges them.
#
# The approximate positions come from the same kernel filter code as the emission
# fast path (`_approx_node_point` / `_exact_node_dir(False(), …)`), minus the
# certification — design §2.4/§2.6, shared code.

@inline function _uf_find(parent::Vector{Int32}, i::Int32)
    while parent[i] != i
        parent[i] = parent[parent[i]]      # path halving
        i = parent[i]
    end
    return i
end

@inline function _uf_union!(parent::Vector{Int32}, i::Int32, j::Int32)
    ri = _uf_find(parent, i); rj = _uf_find(parent, j)
    ri != rj && (parent[ri] = rj)
    return nothing
end

function _merge_coincident_nodes!(m::Manifold, table::NodeTable{P}, seg_nodes; exact) where {P}
    n = num_nodes(table)
    if n < 2
        #-- no cross-key coincidence possible, but one node is routinely reported
        #-- by several candidate pairs (S1; e.g. an a–b–a spike), so the interior
        #-- lists must still be deduped before ordering/splitting (invariant 3)
        for (seg, nids) in seg_nodes
            seg_nodes[seg] = unique(nids)
        end
        return nothing
    end
    parent = collect(Int32(1):Int32(n))
    #-- tier-2 confirmation sweep (manifold-specific proximity geometry)
    nmerges = _coincidence_sweep!(m, table, parent; exact)
    if nmerges == 0
        #-- the overwhelmingly common case (S1: zero coincidences on real data):
        #-- no ids collapse, so skip the O(n) key-rehash compaction entirely.
        #-- Only the interior lists need deduping (a node reported by several pairs);
        #-- the coordinate cache is sized later by `_ensure_coord_cache!`.
        for (seg, nids) in seg_nodes
            seg_nodes[seg] = unique(nids)
        end
        return nothing
    end

    #-- compact: provisional id -> final id, representative = first-seen member
    remap = Vector{Int32}(undef, n)
    final_keys = NodeKey{P}[]
    root_final = Dict{Int32, Int32}()
    for i in 1:n
        r = _uf_find(parent, Int32(i))
        fid = get(root_final, r, Int32(0))
        if fid == 0
            push!(final_keys, table.keys[i])
            fid = Int32(length(final_keys))
            root_final[r] = fid
        end
        remap[i] = fid
    end

    #-- re-point the interner so later endpoint interning resolves to final ids,
    #-- then swap in the compacted key table and a fresh coordinate cache
    for (k, oldid) in table.ids
        table.ids[k] = remap[oldid]
    end
    table.keys = final_keys
    #-- coords/realized are sized once after splitting (`_ensure_coord_cache!`)

    #-- rewrite the interior node lists through the remap, deduping ids a merge
    #-- collapsed together (order is re-established in `split.jl`)
    for (seg, nids) in seg_nodes
        seen = Set{Int32}()
        out = Int32[]
        for id in nids
            fid = remap[id]
            (fid in seen) && continue
            push!(seen, fid); push!(out, fid)
        end
        seg_nodes[seg] = out
    end
    return nothing
end

# Planar proximity sweep: a single x-sort, forward scan while the x-gap is within
# the summed proximity radii, exact confirm, union. A true coincidence has
# near-identical float positions (≪ the crossing radius), so the sweep never
# misses one; over-broad radii only cost extra exact confirms.
function _coincidence_sweep!(m::Planar, table::NodeTable{P}, parent; exact) where {P}
    n = num_nodes(table)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    rad = Vector{Float64}(undef, n)
    for i in 1:n
        k = table.keys[i]
        x, y, _ = _approx_node_point(k)
        xs[i] = x; ys[i] = y
        #-- crossings can land off their float approximation; vertices are exact.
        #-- 1e-8·|coord| is a generous proximity gate (the exact test confirms)
        rad[i] = k.is_crossing ? 1e-8 * max(1.0, abs(x), abs(y)) : 0.0
    end
    order = sortperm(xs)
    nmerges = 0
    @inbounds for a in 1:n
        i = order[a]
        for b in (a + 1):n
            j = order[b]
            rr = rad[i] + rad[j]
            xs[j] - xs[i] > rr && break
            abs(ys[i] - ys[j]) > rr && continue
            _uf_find(parent, Int32(i)) == _uf_find(parent, Int32(j)) && continue
            if rk_nodes_coincide(m, table.keys[i], table.keys[j]; exact)
                _uf_union!(parent, Int32(i), Int32(j)); nmerges += 1
            end
        end
    end
    return nmerges
end

# Spherical proximity sweep: same shape over the (float, normalized) crossing
# directions, gated by chordal proximity along the first coordinate.
function _coincidence_sweep!(m::Spherical, table::NodeTable{P}, parent; exact) where {P}
    n = num_nodes(table)
    dirs = Vector{NTuple{3, Float64}}(undef, n)
    for i in 1:n
        d = _exact_node_dir(False(), table.keys[i])
        s = sqrt(d[1]^2 + d[2]^2 + d[3]^2)
        dirs[i] = (d[1] / s, d[2] / s, d[3] / s)
    end
    xs = [d[1] for d in dirs]
    #-- 1e-11 chord ≈ 6e-4 m on the unit sphere: generous vs the ≤1.4e-14° float
    #-- direction error, so no true coincidence is missed
    rr = 1e-11
    order = sortperm(xs)
    nmerges = 0
    @inbounds for a in 1:n
        i = order[a]
        for b in (a + 1):n
            j = order[b]
            xs[j] - xs[i] > rr && break
            di = dirs[i]; dj = dirs[j]
            dotp = di[1] * dj[1] + di[2] * dj[2] + di[3] * dj[3]
            cx = di[2] * dj[3] - di[3] * dj[2]
            cy = di[3] * dj[1] - di[1] * dj[3]
            cz = di[1] * dj[2] - di[2] * dj[1]
            (dotp > 0 && sqrt(cx^2 + cy^2 + cz^2) <= rr) || continue
            _uf_find(parent, Int32(i)) == _uf_find(parent, Int32(j)) && continue
            if rk_nodes_coincide(m, table.keys[i], table.keys[j]; exact)
                _uf_union!(parent, Int32(i), Int32(j)); nmerges += 1
            end
        end
    end
    return nmerges
end
