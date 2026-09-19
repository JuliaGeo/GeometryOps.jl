# NOTE: This functionality is experimental and may change at any time.

# # Winding-number overlay

#=
## What is a winding-number overlay?

The four [`OverlayNG`](@ref) operations are *binary*: two inputs, A and B, and a
result selected from each face's (location in A, location in B) label pair. A
winding-number overlay is the *N-ary* counterpart. It takes any number of closed
rings — which need not be valid, need not be disjoint, and may self-intersect
freely — nodes all of them together in one pass, and assigns every face of the
resulting arrangement an integer **winding number**: the number of times the
input linework wraps around that face, counting a ring as `+1` on its left and
`-1` on its right.

One predicate on that integer then selects the result:

| `keep` | Operation |
|:--|:--|
| `w -> w > 0` | union of the regions bounded by the rings (an N-ary "unary union") |
| `w -> w >= 1` | the same, and what [`buffer`](@ref) uses |
| `w -> w >= 2` | intersection |
| `isodd` | symmetric difference |

This is the selection rule of Chen & McMains, *Polygon offsetting by computing
winding numbers* (ASME IDETC/CIE 2005), which is why the native planar
[`ChenMcMains`](@ref) buffer is built on it: the raw offset curve of a polygon
is a single self-intersecting ring, and the buffer is exactly its
winding-number-`>= 1` region.

## Why a binary union cannot do this

`union(a, b)` contracts on valid input, and a raw offset curve is a
self-intersecting ring — so it is never a legal operand. Cascading pairwise
unions over a decomposition into valid pieces does produce the right region, at
3-20x the cost, with invalid intermediate results on real input and outright
failures on some of it. The winding overlay takes the self-intersecting ring
directly, which is the whole point of it.

## Implementation

Everything except ring ingest and the seeding step is shared with the binary
engine, and is therefore already manifold-generic:

1. **Ingest** — normalize each ring (drop repeated points, reject non-finite
   coordinates, drop rings with fewer than three distinct vertices, close),
   order the rings along a Hilbert curve so the natural-order segment index
   prunes, and wrap each as a `DIM_A` [`RelateSegmentString`](@ref).
2. **Noding** — the stock `_noded_arrangement`, with `self_node_all = true`:
   the rings are area-contributing but *not* valid areas, so they need the
   all-pairs self-noding pass that `collect.jl` otherwise reserves for linear
   input.
3. **Edge merge and graph** — the stock `_overlay_graph_with_merged_edges`. Each
   ring is given a synthetic `EdgeSourceInfo` with `depth_delta = -1`, so a
   merged edge's summed `a_depth_delta` is exactly `w(right) - w(left)` over
   every ring that shares it. Coincident and opposed rings therefore cancel for
   free.
4. **Labelling** — face cycles from `_face_successor`, one seed per connected
   component, and a flood that propagates `W(left of e) = W(right of e) - d(e)`.
   This is the only manifold-specific step (`_winding_at`); only `Planar()` is
   implemented today, and the shape of the spherical version is recorded below.
5. **Marking and build** — a half-edge is on the result boundary iff
   `keep(W(right)) && !keep(W(left))`, after which `_check_result_area_balance`
   and the stock `_build_polygons` run unchanged.

### Seeding

The binary engine labels faces by locating them against A and B. There is no A
and B here, so each connected component of the arrangement needs one face whose
winding number is known absolutely:

- **Seed node.** The component's lexicographically lowest `(y, x)` node. It is
  always a plain input vertex, never a crossing: a point strictly inside a
  segment has that segment's lower endpoint in the same component, so it cannot
  be the lowest. The first edge of that node's CCW-ordered star has the
  component's outer face on its right.
- **Seed value.** The winding number of that outer face, from a downward
  vertical ray at the seed point, visiting only rings whose bounding box
  contains it. Segments *through* the point are skipped — they are exactly the
  component's own segments at the seed — so the count is well defined.
- **Single component.** Zero, with no ray cast at all.

Components are derived from the arrangement's geometry, not from ring identity,
so duplicated and cancelling rings need no special handling.

### The spherical case (not implemented)

A sphere has no unbounded face and no lowest vertex, so step 4's seeding is the
one part that does not carry over. The design is: pick a reference point `R`
with a known winding number (for a positive buffer distance, any point farther
than the distance from the input has `W_R = 0`; a pole, or the antipode of the
input's bounding-cap centre), count *signed* crossings of the arc `R -> v` for
one vertex `v` per component, and locate the wedge at `v` containing the
direction to `R`. `count_arc_segment!` in `relateng/indexed_point_in_area.jl`
already has the S2-style vertex handling, but counts parity rather than sign.
Everything else — noding, the coincidence sweep, edge merge, star ordering,
cycles, the flood, marking, and `_build_polygons` — is manifold-dispatched
already.
=#

# The winding predicate `buffer` uses, as a named function rather than a closure,
# so the engine specializes on one type across every call.
_winding_at_least_one(w::Integer) = w >= 1

_winding_supported(::Planar) = true
_winding_supported(::Manifold) = false

#=
Normalize one ring to a closed, repeat-free coordinate vector, or `nothing` when
it carries fewer than three distinct vertices (such a ring bounds no area and
contributes nothing to any winding number). Non-finite coordinates are rejected
rather than propagated: they would silently poison every predicate downstream.
=#
function _winding_normalize_ring(::Type{P}, ring) where {P <: Tuple{Float64, Float64}}
    out = P[]
    sizehint!(out, length(ring) + 1)
    for p in ring
        q = (Float64(p[1]), Float64(p[2]))
        (isfinite(q[1]) && isfinite(q[2])) ||
            throw(ArgumentError("winding overlay: non-finite ring coordinate $q"))
        (isempty(out) || out[end] != q) && push!(out, q)
    end
    while length(out) > 1 && out[end] == out[1]
        pop!(out)
    end
    length(out) < 3 && return nothing
    #-- three or more points but only two distinct ones is a degenerate a-b-a spike
    a = out[1]; b = out[2]
    any(p -> p != a && p != b, out) || return nothing
    push!(out, out[1])
    return out
end

# The ring bounding boxes, in the same order as the rings. Planar only: the
# seeding ray and the Hilbert order are both planar notions.
_winding_ring_extents(::Planar, rings::Vector{Vector{P}}) where {P} =
    [_winding_ring_extent(r) for r in rings]

function _winding_ring_extent(pts::Vector{P}) where {P}
    xlo = xhi = pts[1][1]; ylo = yhi = pts[1][2]
    @inbounds for (x, y) in pts
        xlo = min(xlo, x); xhi = max(xhi, x)
        ylo = min(ylo, y); yhi = max(yhi, y)
    end
    return Extents.Extent(X = (xlo, xhi), Y = (ylo, yhi))
end

#=
Rings in Hilbert order of their bounding-box centres. The segment index the
noder builds is a natural-order (`Unsorted`) `RTree`, so it prunes well only
when consecutive segments are spatially close — which holds inside a ring and
fails across a spatially incoherent list of rings. Measured on 10 000 scattered
octagons, this reordering took the arrangement from 3306 ms to 87 ms.

Reuses the `HPR` loader's key function, so the buffer and the R-trees quantize
space the same way.
=#
function _winding_spatial_order(extents::Vector{<:Extents.Extent})
    length(extents) <= 1 && return collect(eachindex(extents))
    return sortperm(FlexibleRTrees._hilbert_keys(extents))
end

"""
    _winding_overlay(m::Manifold, ::Type{T}, rings; keep, exact = True())

The N-ary winding-number overlay: node every ring in `rings` together, label
each face of the arrangement with its winding number, and return the faces
satisfying `keep` as a `Vector` of polygons over output point type `T`.

`rings` is a vector of closed coordinate vectors in the manifold's kernel point
type. They may self-intersect, coincide, and repeat; each contributes `+1` to
the winding number on its left. `keep` receives an `Int` and should be a named
function or a callable struct, not an anonymous closure built per call, so the
engine specializes once.

Only `Planar()` is implemented; see the file header for the spherical design.
"""
function _winding_overlay(m::Manifold, ::Type{T}, rings::AbstractVector;
        keep::K = _winding_at_least_one, exact = True()) where {T, K}
    _winding_supported(m) || throw(ArgumentError(
        "the winding-number overlay is implemented on the `Planar()` manifold " *
        "only; got $(typeof(m)). See `winding_overlay.jl` for the spherical design."))
    P = _kernel_point_type(m)
    #-- 1. ingest
    norm = Vector{P}[]
    for r in rings
        pts = _winding_normalize_ring(P, r)
        pts === nothing || push!(norm, pts)
    end
    isempty(norm) && return _result_poly_type(T)[]
    extents = _winding_ring_extents(m, norm)
    order = _winding_spatial_order(extents)
    ss = [_rss_create_ring(norm[i], true, j, 1, nothing, nothing) for (j, i) in enumerate(order)]
    #-- the per-ring index used for seeding; `STR` because it is queried, not traversed
    ring_tree = RTree(STR(), collect(eachindex(ss)); extents = extents[order])
    #-- built here and handed to the noder, which would otherwise build it twice
    #-- (once for the A x B pass it never runs, once for self-noding). Monotone
    #-- chains where the manifold allows them (`chains.jl`), segments otherwise.
    seg_tree = _noding_index(m, ss)

    #-- 2. noding. `self_node_all`: these rings are area-contributing but not
    #-- valid areas, so the vertex-scope self-noding pass is not sufficient.
    #-- Splitting is `run_split.jl`'s, not `split.jl`'s: one graph edge per
    #-- polyline run between crossings, not one per curve segment.
    arr, vids, voff = _winding_arrangement(m, T, ss, seg_tree; exact)

    #-- 3. edge merge + graph, keeping the summed integer deltas.
    #-- `depth_delta = -1` per ring: `EdgeSourceInfo`'s sign convention is
    #-- `w(right) - w(left)` (see `_location_left`), and a ring adds +1 on its LEFT.
    sources = fill(EdgeSourceInfo(Int8(0), DIM_A, false, Int8(-1)), length(ss))
    merged = _merge_noded_runs(arr, sources)
    g = _graph_from_merged(m, arr, merged; exact)
    runs = WindingRuns(merged, vids, voff)

    #-- 4. winding labels, 5. rings -> polygons (both stock from here)
    marked = _label_winding!(m, g, merged, ring_tree, ss, keep; exact)
    _check_result_area_balance(g)
    return _build_polygons(m, g, marked; exact, runs)
end

#=
The winding overlay's noding pass: `_noded_arrangement`'s stages 1 and 3 verbatim
(collect, then coincidence merging), with `run_split.jl` in place of `split.jl`
for stages 2 and 4. Returns the arrangement plus the per-position vertex-id table
the ring builder needs to emit a run's interior vertices.
=#
function _winding_arrangement(m::Manifold, ::Type{T},
        ss::Vector{RelateSegmentString{P}}, seg_tree; exact) where {P, T}
    segstrings = copy(ss)
    table = NodeTable{P, T}()
    seg_nodes = NTuple{3, Int32}[]
    _collect_crossings!(m, table, seg_nodes, ss, RelateSegmentString{P}[], Int32(length(ss));
                        exact, tree_a = seg_tree, tree_b = nothing, clip_a = nothing,
                        clip_b = nothing, self_node_all = true)
    _merge_coincident_nodes!(m, table, seg_nodes; exact)
    edges, vids, voff = _split_runs!(m, table, seg_nodes, segstrings; exact)
    _ensure_coord_cache!(table)
    arr = NodedArrangement{P, T}(segstrings, table, seg_nodes, edges, BitVector())
    return arr, vids, voff
end

#=
Winding number at `p`, counted along a downward vertical ray, over the rings
whose bounding box contains `p` (every other ring contributes 0). `p` must lie
on no counted segment; segments through it (`orient == 0`) are skipped, which is
exactly the caller's own component at its seed vertex.

The x-rule is half-open so a vertex shared by two segments is counted once: a
rightward segment with `p` on its left adds +1, a leftward segment with `p` on
its right adds -1. Orientation is the exact kernel predicate, so the count is
exact even for a ray grazing a near-horizontal segment.
=#
function _winding_at(::Planar, ring_tree, ss, p::P) where {P}
    px, py = p
    w = 0
    pred = ext -> (ext.X[1] <= px <= ext.X[2]) & (ext.Y[1] <= py <= ext.Y[2])
    SpatialTreeInterface.depth_first_search(pred, ring_tree) do i
        w += _winding_ray_count(ss[ring_tree.data[i]].pts, px, py, p)
        return nothing
    end
    return w
end

function _winding_ray_count(pts::Vector{P}, px, py, p) where {P}
    w = 0
    @inbounds for k in 1:(length(pts) - 1)
        a = pts[k]; b = pts[k + 1]
        (a[2] > py && b[2] > py) && continue
        if a[1] <= px < b[1]            # rightward, half-open at the right end
            rk_orient(Planar(), a, b, p; exact = True()) > 0 && (w += 1)
        elseif b[1] <= px < a[1]        # leftward
            rk_orient(Planar(), a, b, p; exact = True()) < 0 && (w -= 1)
        end
    end
    return w
end

#=
Label every face with its winding number and mark the half-edges on the result
boundary, returning the marked indices for `_build_polygons`.

The per-half-edge delta comes from `merged` through the pairing the graph builder
documents and asserts: half-edges `2i-1` (forward) and `2i` (reverse) are the two
halves of `merged[i]`, and the reverse half sees the negated delta.
=#
function _label_winding!(m::Manifold, g::OverlayGraph, merged::Vector{MergeEdge},
        ring_tree, ss, keep::K; exact) where {K}
    E = g.edges
    nE = length(E)
    keys = g.arr.nodes.keys
    @assert nE == 2 * length(merged) """
        winding overlay: the graph holds $nE half-edges for $(length(merged)) merged \
        edges; the forward/reverse pairing this labeller reads is broken"""
    delta = Vector{Int32}(undef, nE)
    @inbounds for i in eachindex(merged)
        delta[2i - 1] = merged[i].a_depth_delta
        delta[2i] = -merged[i].a_depth_delta
    end

    #-- face cycles: cyc[e] identifies the face on the RIGHT of half-edge e
    cyc = zeros(Int32, nE)
    ncyc = Int32(0)
    for i in 1:nE
        cyc[i] != 0 && continue
        ncyc += Int32(1)
        e = Int32(i)
        while true
            cyc[e] = ncyc
            e = _face_successor(E, e)
            e == i && break
        end
    end

    #-- connected components of the arrangement (union-find over merged edges)
    nn = length(keys)
    uf = collect(Int32(1):Int32(nn))
    for me in merged
        _uf_union!(uf, me.node_lo, me.node_hi)
    end

    #-- one seed per component: its lexicographically lowest (y, x) NODE, which is
    #-- always a plain input vertex (see the header), hence never a crossing key
    lowest = Dict{Int32, Int32}()
    for nid in Int32(1):Int32(nn)
        g.node_edges[nid] == 0 && continue
        k = keys[nid]
        k.is_crossing && continue
        r = _uf_find(uf, nid)
        cur = get(lowest, r, Int32(0))
        if cur == 0 || (k.pt[2], k.pt[1]) < (keys[cur].pt[2], keys[cur].pt[1])
            lowest[r] = nid
        end
    end

    W = fill(typemin(Int32), ncyc)
    todo = Int32[]
    single = length(lowest) == 1
    for (_, v) in lowest
        #-- the smallest-angle out-edge of the star has the component's outer
        #-- face on its right
        c0 = cyc[g.node_edges[v]]
        W[c0] = single ? Int32(0) : Int32(_winding_at(m, ring_tree, ss, keys[v].pt))
        push!(todo, c0)
    end

    #-- one representative half-edge per cycle, then flood across `sym` pairs
    cyc_start = zeros(Int32, ncyc)
    for i in Int32(1):Int32(nE)
        cyc_start[cyc[i]] == 0 && (cyc_start[cyc[i]] = i)
    end
    while !isempty(todo)
        cy = pop!(todo)
        e0 = cyc_start[cy]
        e = e0
        while true
            cs = cyc[he_sym(E, e)]
            wl = W[cy] - delta[e]
            if W[cs] == typemin(Int32)
                W[cs] = wl
                push!(todo, cs)
            elseif W[cs] != wl
                throw(_OverlayTopologyError(
                    "winding overlay: inconsistent winding at half-edge $e — face " *
                    "$cs is $(W[cs]) from one side and $wl from another"))
            end
            e = _face_successor(E, e)
            e == e0 && break
        end
    end

    marked = Int32[]
    for e in Int32(1):Int32(nE)
        W[cyc[e]] == typemin(Int32) && throw(_OverlayTopologyError(
            "winding overlay: face $(cyc[e]) was never labelled"))
        if keep(Int(W[cyc[e]])) && !keep(Int(W[cyc[he_sym(E, e)]]))
            oe_mark_in_result_area!(E, e)
            push!(marked, e)
        end
    end
    return marked
end
