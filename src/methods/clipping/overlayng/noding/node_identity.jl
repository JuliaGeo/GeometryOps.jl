# NOTE: This functionality is experimental and may change at any time.

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

# Sort the flat `(string, segment, node)` record and drop exact duplicates: one
# node is routinely reported by several candidate pairs (S1; e.g. an a–b–a
# spike), and `split.jl` requires each segment's interior ids distinct
# (invariant 3) and grouped by segment. `split.jl` re-establishes the geometric
# order within a group, so sorting by node id costs nothing.
function _sort_seg_nodes!(v::Vector{NTuple{3, Int32}})
    sort!(v)
    isempty(v) && return v
    w = 1
    @inbounds for r in 2:length(v)
        if v[r] != v[w]
            w += 1
            v[w] = v[r]
        end
    end
    resize!(v, w)
    return v
end

function _merge_coincident_nodes!(m::Manifold, table::NodeTable{P},
        seg_nodes::Vector{NTuple{3, Int32}}; exact) where {P}
    _sort_seg_nodes!(seg_nodes)
    n = num_nodes(table)
    #-- fewer than two nodes: no cross-key coincidence is possible
    n < 2 && return nothing
    parent = collect(Int32(1):Int32(n))
    #-- tier-2 confirmation sweep (manifold-specific proximity geometry)
    nmerges = _coincidence_sweep!(m, table, parent; exact)
    #-- the overwhelmingly common case (S1: zero coincidences on real data): no
    #-- ids collapse, so skip the O(n) key-rehash compaction entirely. The
    #-- coordinate cache is sized later by `_ensure_coord_cache!`.
    nmerges == 0 && return nothing

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

    #-- rewrite the interior node record through the remap, re-deduping ids a
    #-- merge collapsed together (order is re-established in `split.jl`)
    @inbounds for i in eachindex(seg_nodes)
        (g, k, id) = seg_nodes[i]
        seg_nodes[i] = (g, k, remap[id])
    end
    _sort_seg_nodes!(seg_nodes)
    return nothing
end

# Planar proximity sweep: a single x-sort, then a cascade of interval groupings —
# proximity radii, then certified radii — bottoming out in an exact-point
# bucketing. Every node is positioned to within max(1e-9·scale, half an ulp) — a
# crossing whose certified float error exceeds 1e-9·scale is re-positioned from
# its exact rational point — so the 1e-8·scale gate strictly exceeds any
# coincident pair's summed error and no merge can be missed; over-broad radii
# only cost extra exact confirms. This gate is sound only because
# `_approx_node_point`'s radius is conservative.
#=
## Why the pairwise sweep is not enough (the degenerate-node fix)

A forward pairwise scan is O(#pairs the x-break survives), and the break cannot
separate nodes that share an x. `k` concurrent lines produce `C(k, 2)` crossing
keys at one point, so the scan degenerates to `C(k,2)²/2` pairs, and — because
those crossings almost never coincide *exactly* in Float64 — every one of them
falls through to `rk_nodes_coincide`, which builds two uncached
`Rational{BigInt}` crossing points (~20 µs). Eroding a regular `n`-gon by exactly
its inradius makes all `n` offset lines concurrent; the sweep then cost Θ(n⁴)
exact constructions, all of them returning `false`: 2.3 s at `n = 32`, 281 s at
`n = 100`.

So the pairwise test is replaced by a grouping. `rk_nodes_coincide` *is* equality
of `_exact_node_point`, and equality classes are found by hashing, not by
comparing pairs:

1. connected components of the proximity intervals `[x ± rad]`, then of `[y ± rad]`
   within each — the transitive closure of today's pair test, so no candidate is
   lost;
2. the same two splits again on the *certified* radius `err` that
   `_approx_node_point` reports, which proves most nearby-but-distinct nodes apart
   without any exact arithmetic;
3. inside each surviving group, one `_exact_node_point` per member into a `Dict`,
   unioning members that land in the same bucket.

Exact constructions go from Θ(#pairs) to Θ(#nodes in a non-singleton group), and
the group walk from quadratic to `O(n log n)`. Singletons — every node on real
data, where the S1 census found zero coincidences — cost nothing but the sort.
=#

# Connected components of the intervals `[v[i] - r[i], v[i] + r[i]]` over `idx`,
# which must already be sorted by `v`. Calls `f(a, b)` with each component's
# inclusive slice of `idx`. Interval overlap is exactly the pair test
# `|v[i] - v[j]| <= r[i] + r[j]`, so a component is its transitive closure.
function _interval_runs(f::F, idx, v::Vector{Float64}, r::Vector{Float64}) where {F}
    n = length(idx)
    a = 1
    @inbounds while a <= n
        b = a
        reach = v[idx[a]] + r[idx[a]]
        while b < n && v[idx[b + 1]] - r[idx[b + 1]] <= reach
            b += 1
            reach = max(reach, v[idx[b]] + r[idx[b]])
        end
        f(a, b)
        a = b + 1
    end
    return nothing
end

# Sort `mem` (a candidate group) on `v`, split it into components at radius `r`,
# and hand each component of two or more members to `next`.
function _split_group!(mem::Vector{Int32}, v::Vector{Float64}, r::Vector{Float64},
        next::F) where {F}
    length(mem) < 2 && return 0
    sort!(mem; by = i -> @inbounds v[i])
    total = 0
    _interval_runs(mem, v, r) do a, b
        b > a && (total += next(mem[a:b]))
        return nothing
    end
    return total
end

function _coincidence_sweep!(m::Planar, table::NodeTable{P}, parent; exact) where {P}
    n = num_nodes(table)
    xs = Vector{Float64}(undef, n)
    ys = Vector{Float64}(undef, n)
    rad = Vector{Float64}(undef, n)
    #-- the certified radius of the float estimate, used to prove groups apart
    cert = Vector{Float64}(undef, n)
    EP = Tuple{Rational{BigInt}, Rational{BigInt}}
    epts = Vector{EP}(undef, n)
    have = falses(n)
    for i in 1:n
        k = table.keys[i]
        x, y, err = _approx_node_point(k)
        #-- negated so a NaN solve (float determinant rounded to 0) also escalates
        if k.is_crossing && !(err <= 1e-9 * max(1.0, abs(x), abs(y)))
            p = _exact_node_point(k)
            epts[i] = p; have[i] = true
            x = Float64(p[1]); y = Float64(p[2]); err = 0.0
        end
        xs[i] = x; ys[i] = y; cert[i] = err
        #-- crossings can land off their float approximation; vertices are exact.
        #-- 1e-8·|coord| is a generous proximity gate (the exact test confirms)
        rad[i] = k.is_crossing ? 1e-8 * max(1.0, abs(x), abs(y)) : 0.0
    end

    #-- the bottom of the cascade: one exact point per member, bucketed
    buckets = Dict{EP, Int32}()
    nmerges = 0
    exact_group! = function (mem::Vector{Int32})
        length(mem) < 2 && return 0
        empty!(buckets)
        merged = 0
        @inbounds for i in mem
            have[i] || (epts[i] = _exact_node_point(table.keys[i]); have[i] = true)
            rep = get(buckets, epts[i], Int32(0))
            if rep == Int32(0)
                buckets[epts[i]] = i
            elseif _uf_find(parent, rep) != _uf_find(parent, i)
                _uf_union!(parent, rep, i); merged += 1
            end
        end
        return merged
    end
    #-- certified split (x then y) before any exact arithmetic is spent
    cert_y! = mem -> _split_group!(mem, ys, cert, exact_group!)
    cert_x! = mem -> _split_group!(mem, xs, cert, cert_y!)
    prox_y! = mem -> _split_group!(mem, ys, rad, cert_x!)

    order = Int32[Int32(i) for i in sortperm(xs)]
    _interval_runs(order, xs, rad) do a, b
        b > a && (nmerges += prox_y!(order[a:b]))
        return nothing
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
