# # OverlayNG noding substrate — the `NodedArrangement`
#
# Phase 1 of the OverlayNG port (design doc `2026-07-16-overlayng-noding-substrate.md`).
# Nodes two input geometries into a shared, exactly-noded edge arrangement with
# **symbolic** node identity: node identity, ordering, and coincidence are decided
# by exact kernel predicates over the input vertices and the symbolic crossing
# keys (`NodeKey`); Float64 appears only as certified filters and at
# emission. There is no snapping, no tolerance in any decision (design §0).
#
# The pipeline is four stages over the reused RelateNG substrate:
#   1. collect   (`collect.jl`)      — candidate enumeration + exact classification
#   2. identity  (`node_identity.jl`) — two-tier node grouping
#   3. order+split (`split.jl`)       — along-segment ordering, noded-edge emission
#   4. emit      (`emit.jl`)          — the sole lossy step, coordinate realization
#
# Everything here is internal to GeometryOps — nothing is exported.

# A noded sub-segment: a piece of one parent segment between two nodes. It
# carries no geometry (design §2.1 invariant 5) — its shape is a lookup into its
# parent segment, and all source metadata (owner, ring id, dimension) is reached
# through `string_idx` into the arrangement's `segstrings`, so nothing can
# desynchronize.
struct NodedEdge
    string_idx :: Int32    # index into `NodedArrangement.segstrings`
    seg_idx    :: Int32    # segment within that parent string (`pts[seg_idx:seg_idx+1]`)
    node_lo    :: Int32    # node id at the sub-segment start (in the parent's traversal order)
    node_hi    :: Int32    # node id at the sub-segment end
end

#=
The node identity table (design §2.4). `ids` is the tier-1 egal interner
(`NodeKey` bit-equality is canonical — kernel points normalize signed zeros),
mapping every known key to its node id; after tier-2 merging it maps every
provisional key to its *final* id, so keys interned later (segment endpoints in
`split.jl`) resolve to the merged id. `keys[id]` is the group representative.
`coords`/`realized` memoize emitted output coordinates (design §2.6), realized
lazily by `node_point` and grown as endpoint nodes are interned.

Mutable so tier-2 merging can compact `keys` and re-point `ids` in place.
=#
mutable struct NodeTable{P}
    ids      :: Dict{NodeKey{P}, Int32}
    keys     :: Vector{NodeKey{P}}
    coords   :: Vector{Tuple{Float64, Float64}}
    realized :: Vector{Bool}
end

NodeTable{P}() where {P} =
    NodeTable{P}(Dict{NodeKey{P}, Int32}(), NodeKey{P}[], Tuple{Float64, Float64}[], Bool[])

num_nodes(t::NodeTable) = length(t.keys)

# Intern a key, returning its node id (tier-1 egal merge). The output-coordinate
# cache is an emission concern (`_ensure_coord_cache!`), not grown here.
function _intern_node!(t::NodeTable{P}, key::NodeKey{P}) where {P}
    id = get(t.ids, key, Int32(0))
    id != 0 && return id
    push!(t.keys, key)
    id = Int32(length(t.keys))
    t.ids[key] = id
    return id
end

# Size the (lazily-realized) output-coordinate cache to the final node count.
# Called once after splitting, so noding itself never touches it.
function _ensure_coord_cache!(t::NodeTable)
    n = length(t.keys)
    t.coords = Vector{Tuple{Float64, Float64}}(undef, n)
    t.realized = fill(false, n)
    return nothing
end

"""
    NodedArrangement{P}

The exactly-noded arrangement of two input geometries (design §2.1). `P` is the
manifold's kernel point type — exactly two instantiations,
`Tuple{Float64,Float64}` (planar) and `UnitSphericalPoint{Float64}` (spherical) —
so the engine is type-erased over the input geometry types.

Fields:
- `segstrings`: the ingested inputs as `RelateSegmentString`s (A side
  first, then B side); `NodedEdge.string_idx` indexes here.
- `nodes`: the symbolic node table (`NodeTable`).
- `seg_nodes`: per-parent-segment ordered interior node-id lists, keyed by
  `(string_idx, seg_idx)`; absent for unsplit segments.
- `edges`: every noded sub-segment of every parent segment.

Construct with `NodedArrangement(m, a, b)` (raw geometries) or
`NodedArrangement(m, ssa, ssb)` (pre-ingested segment strings).
"""
struct NodedArrangement{P}
    segstrings :: Vector{RelateSegmentString{P}}
    nodes      :: NodeTable{P}
    seg_nodes  :: Dict{Tuple{Int32, Int32}, Vector{Int32}}
    edges      :: Vector{NodedEdge}
end

num_nodes(arr::NodedArrangement) = num_nodes(arr.nodes)
num_edges(arr::NodedArrangement) = length(arr.edges)

# Whether a segment string is a polygon hole. Derived (not stored — design §2.2):
# shells are `ring_id == 0`, holes `ring_id >= 1`, lines `ring_id == -1`
# (`_extract_ring_to_segment_string!` in relate_geometry.jl).
_ss_is_hole(ss::RelateSegmentString) = ss.ring_id > 0

# ## Ingest (design §2.2)
#
# Reuse `RelateGeometry` / `extract_segment_strings` unchanged: kernel-point
# conversion, repeated-point removal, and ring orientation all happen here, once.
# Kept separate from arrangement construction so a future prepared overlay can
# convert once and re-arrange many times (S1: extraction dominates sparse pairs).
_overlay_segstrings(m::Manifold, geom, is_a::Bool; exact = True()) =
    extract_segment_strings(RelateGeometry(m, geom; exact), is_a, nothing)

# ## Construction

# From raw geometries: ingest each side, then arrange.
function NodedArrangement(m::Manifold, a, b; exact = True(), tree_a = nothing, tree_b = nothing)
    ssa = _overlay_segstrings(m, a, true; exact)
    ssb = _overlay_segstrings(m, b, false; exact)
    return NodedArrangement(m, ssa, ssb; exact, tree_a, tree_b)
end

# From pre-ingested segment strings. `tree_a`/`tree_b` accept caller-supplied
# prebuilt segment indices (a `PreparedRelate` carries exactly `_relate_edge_index`
# output) — an optional argument only, no new prepare type (design §2.3).
function NodedArrangement(m::Manifold,
        ssa::AbstractVector{RelateSegmentString{P}},
        ssb::AbstractVector{RelateSegmentString{P}};
        exact = True(), tree_a = nothing, tree_b = nothing) where {P}
    na = length(ssa)
    segstrings = Vector{RelateSegmentString{P}}(undef, na + length(ssb))
    @inbounds for i in 1:na
        segstrings[i] = ssa[i]
    end
    @inbounds for i in eachindex(ssb)
        segstrings[na + i] = ssb[i]
    end

    table = NodeTable{P}()
    seg_nodes = Dict{Tuple{Int32, Int32}, Vector{Int32}}()
    # stage 1
    _collect_crossings!(m, table, seg_nodes, ssa, ssb, Int32(na); exact, tree_a, tree_b)
    # stage 3
    _merge_coincident_nodes!(m, table, seg_nodes; exact)
    # stages 2 + 4
    edges = _split_edges!(m, table, seg_nodes, segstrings; exact)
    _ensure_coord_cache!(table)
    return NodedArrangement{P}(segstrings, table, seg_nodes, edges)
end
