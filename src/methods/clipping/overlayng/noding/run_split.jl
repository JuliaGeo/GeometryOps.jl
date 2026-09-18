# NOTE: This functionality is experimental and may change at any time.

# # Stages 2 + 4 for the winding overlay — split into POLYLINE runs
#
# A graph edge here is a *run*: the polyline a parent string traces between two
# consecutive nodes. This is JTS's shape — a noded substring splits at
# intersection nodes and nowhere else — and it sizes the graph by the crossings
# rather than by the vertices: a crossing-free 20 000-segment offset curve
# becomes three nodes and three edges, and its vertices ride along as a
# coordinate range on the run's parent.
#
# `split.jl` keeps the one-edge-per-segment shape the binary overlay wants, where
# a graph edge IS a segment and every vertex is a node.
#
# ## What a "stop" is
#
# A parent position becomes a graph node (a *stop*) when it is
#
#   1. a node the noder already created — a proper crossing on the segment, or a
#      vertex the noder found incident on foreign linework (`table.ids`);
#   2. a vertex coordinate shared by two parent positions. The noder leaves these
#      to the caller: `_classify_pair!` guards every vertex record with
#      `a0 != b0 && a0 != b1`, so a vertex sitting exactly on another curve's
#      *vertex* is recorded nowhere, and `split.jl` merges them implicitly by
#      interning every endpoint into one table. The coordinate pass below does it
#      explicitly — a vertex seen twice is promoted to a stop; or
#   3. one of three forced positions per closed string.
#
# ## The three forced stops
#
# Every closed string forces a stop at
#
#   * **position 1**, so a run never wraps the closure and `seg_lo <= seg_hi`
#     always holds;
#   * the **`(y, x)`-minimal vertex**, which `_label_winding!` depends on: it
#     seeds each connected component at the component's lexicographically lowest
#     NODE and needs the smallest-angle star edge there to face the outer face.
#     That argument is about the lowest POINT of the component, so the lowest
#     point has to be a graph node — and it is always a vertex of some string,
#     never a crossing (a crossing lies strictly inside a segment, whose lower
#     endpoint precedes it);
#   * a **midpoint**, so a string always carries two distinct stops. One stop
#     would make the whole string a single self-loop edge, whose
#     `node_lo == node_hi` defeats the merge direction test (`_merge` reads
#     `inc.node_lo == base.node_lo`).
#
# Three interned keys per string, where `split.jl` interns one per segment.

#=
Per-string vertex identity and the runs' interior coordinate ranges.

`vids[voff[s] + j]` is the node id of vertex position `j` of string `s`
(`1 <= j <= length(pts)-1`; a closed string's position `length(pts)` IS position
1). `merged[mi]` carries the run's parent and its `seg_idx : seg_hi` segment
span, whose interior vertices are positions `seg_idx+1 : seg_hi`.

Held apart from `NodedArrangement` because only the winding path produces one,
and the binary overlay must keep reading the arrangement it has always read.
=#
struct WindingRuns
    merged :: Vector{MergeEdge}
    vids   :: Vector{Int32}
    voff   :: Vector{Int32}
end

# The interior node ids of merged run `mi`, walked in the half-edge's direction.
@inline function _run_interior_range(r::WindingRuns, mi::Integer)
    me = @inbounds r.merged[mi]
    base = @inbounds r.voff[me.string_idx]
    return (base + me.seg_idx + Int32(1)):(base + me.seg_hi)
end

# ## Splitting

#=
Split every (closed) parent string into polyline runs between consecutive stops.

Returns `(edges, vids, voff)`. `table` gains one `NodeKey` per non-stop vertex —
appended directly to `table.keys` with an arithmetically assigned id and NOT
interned into `table.ids`, because a non-stop vertex is by construction unique:
had its coordinate been reached by any other position, rule 2 above would have
made it a stop. Those ids exist only so the ring builder and the exact
orientation/area predicates keep seeing a node id per emitted coordinate.

Almost every graph node id is `<= num_nodes(table)` as of the forced-stop pass;
the exceptions are plain vertices promoted by rule 2, which are rare and which
the CSR star array covers anyway (it is sized by the full node count, and a
plain vertex simply has an empty slice).
=#
function _split_runs!(m::Manifold, table::NodeTable{P},
        seg_nodes::Vector{NTuple{3, Int32}},
        segstrings::Vector{RelateSegmentString{P}}; exact) where {P}
    nS = length(segstrings)
    #-- 0. per-string base offset into the flat vertex-id array
    voff = Vector{Int32}(undef, nS + 1)
    tot = Int32(0)
    @inbounds for s in 1:nS
        voff[s] = tot                     # position j lives at voff[s] + j
        tot += Int32(length(segstrings[s].pts) - 1)
    end
    voff[nS + 1] = tot

    #-- 1. force the three stops per string, BEFORE the identity pass, so every
    #--    forced coordinate is already in the table when the pass reads it
    for s in 1:nS
        pts = segstrings[s].pts
        n = length(pts) - 1
        _intern_node!(table, vertex_node(pts[1]))
        jmin = 1
        @inbounds for j in 2:n
            p = pts[j]; q = pts[jmin]
            ((p[2], p[1]) < (q[2], q[1])) && (jmin = j)
        end
        jmin != 1 && _intern_node!(table, vertex_node(pts[jmin]))
        #-- a second distinct stop: the midpoint, or the first position that is
        #-- not a repeat of position 1 if the midpoint happens to be one
        jmid = 1 + (n >> 1)
        if jmid <= n && pts[jmid] != pts[1]
            _intern_node!(table, vertex_node(pts[jmid]))
        elseif jmin == 1
            @inbounds for j in 2:n
                if pts[j] != pts[1]
                    _intern_node!(table, vertex_node(pts[j])); break
                end
            end
        end
    end
    ngraph = Int32(num_nodes(table))

    #-- 2. vertex identity, one coordinate hash per position. Seeded from the
    #--    noder's own vertex keys (post-merge `ids` maps every interned key to
    #--    its FINAL id, so a vertex merged into a crossing resolves correctly).
    vdict = Dict{P, Int32}()
    sizehint!(vdict, length(table.ids) + tot)
    for (k, id) in table.ids
        k.is_crossing || (vdict[k.pt] = id)
    end
    #-- every position either resolves to an existing node or creates one, so
    #-- both arrays are sized once here. Growing `table.keys` by doubling instead
    #-- costs ~2x the final byte count in memcpy, and a `NodeKey` is 72 bytes:
    #-- 57 MB of it on the 800k-segment meander.
    nkeys = table.keys
    sizehint!(nkeys, Int(ngraph) + Int(tot))
    is_stop = Vector{Bool}(undef, Int(ngraph) + Int(tot))
    fill!(view(is_stop, 1:Int(ngraph)), true)
    vids = Vector{Int32}(undef, tot)
    for s in 1:nS
        pts = segstrings[s].pts
        base = voff[s]
        @inbounds for j in 1:(length(pts) - 1)
            p = pts[j]
            id = get(vdict, p, Int32(0))
            if id == 0
                push!(nkeys, vertex_node(p))
                id = Int32(length(nkeys))
                is_stop[id] = false
                vdict[p] = id
            elseif !is_stop[id]
                #-- second sighting of a plain vertex coordinate: a junction
                is_stop[id] = true
            end
            vids[base + j] = id
        end
    end

    #-- 3. order each split segment's interior nodes along that segment. The flat
    #--    record is sorted by `(string, segment, node)`, so its `(string,
    #--    segment)` groups are already contiguous: one pass finds each, and only
    #--    a group of two or more reaches the kernel comparator. The ordered ids
    #--    are written back, which is what leaves the record in geometric rather
    #--    than id order within a group — the same state `split.jl` hands the
    #--    binary path, and step 4 below is its only reader.
    nrec = length(seg_nodes)
    ord = Int32[]
    i = 1
    @inbounds while i <= nrec
        (g, k, _) = seg_nodes[i]
        j = i + 1
        while j <= nrec && seg_nodes[j][1] == g && seg_nodes[j][2] == k
            j += 1
        end
        if j - i >= 2
            resize!(ord, j - i)
            for t in i:(j - 1)
                ord[t - i + 1] = seg_nodes[t][3]
            end
            pts = segstrings[g].pts
            _order_along_segment!(m, ord, pts[k], pts[k + 1], table; exact)
            for t in i:(j - 1)
                seg_nodes[t] = (g, k, ord[t - i + 1])
            end
        end
        i = j
    end

    #-- 4. walk each string, emitting one run per consecutive stop pair. One
    #--    cursor walks the flat record in lockstep with the segment loop, so a
    #--    segment carrying no interior node costs one comparison.
    edges = NodedEdge[]
    #-- one run per stop, and a stop is an interior node or a forced position
    sizehint!(edges, nrec + 3 * nS)
    cur = 1
    for s in 1:nS
        s32 = Int32(s)
        @inbounds while cur <= nrec && seg_nodes[cur][1] < s32
            cur += 1
        end
        pts = segstrings[s].pts
        n = Int32(length(pts) - 1)
        base = voff[s]
        #-- position 1 is forced, so the walk starts and ends on a real stop
        start_id = @inbounds vids[base + 1]
        cur_id = start_id
        cur_seg = Int32(1)                       # first segment of the open run
        @inbounds for k in Int32(1):n
            #-- a stop at the vertex that OPENS segment k closes the run at k-1
            if k > 1 && is_stop[vids[base + k]]
                _emit_run!(edges, s32, cur_seg, k - Int32(1), cur_id, vids[base + k])
                cur_id = vids[base + k]; cur_seg = k
            end
            #-- then the nodes strictly inside segment k, in along-segment order
            while cur <= nrec
                (g, kk, nid) = seg_nodes[cur]
                (g != s32 || kk != k) && break
                #-- `nid == cur_id` is NOT skipped here: a string that crosses
                #-- itself meets the same node at two positions, and the piece
                #-- between them is a genuine (loop) run. Only a run with no
                #-- vertex of its own is dropped, by `_emit_run!`.
                _emit_run!(edges, s32, cur_seg, k, cur_id, nid)
                cur_id = nid; cur_seg = k
                cur += 1
            end
        end
        #-- close back onto position 1 (== position n+1 for a closed string)
        _emit_run!(edges, s32, cur_seg, n, cur_id, start_id)
    end
    return edges, vids, voff
end

@inline function _emit_run!(edges::Vector{NodedEdge}, s::Int32, seg_lo::Int32,
        seg_hi::Int32, lo::Int32, hi::Int32)
    #-- a zero-length link (two stops that coincidence merging collapsed onto one
    #-- id, with no vertex between them) carries no geometry, exactly as the
    #-- `a == b` skip in `split.jl`
    (lo == hi && seg_hi < seg_lo + Int32(1)) && return nothing
    push!(edges, NodedEdge(s, seg_lo, lo, hi, seg_hi))
    return nothing
end

# ## Merging coincident runs
#
# `_merge_noded_edges` keys on the unordered node pair, which is unique for a
# straight sub-segment. It is NOT unique for a run: two runs can join the same
# pair of nodes by different paths — the two halves of a crossing-free closed
# curve always do. The key therefore also carries the run's length in vertices
# and its first interior vertex, read in the canonical `(lo -> hi)` direction.
#
# That is enough because a run's interior vertices are never nodes: two runs that
# share an endpoint pair AND a first interior coordinate would have that
# coordinate reached by two positions, which `_split_runs!` rule 2 promotes to a
# stop, which would have split both runs there. And two runs that are coincident
# necessarily have identical vertex sequences, for the same reason — a vertex of
# one lying inside the other is an incidence the noder records.

@inline _run_key_sentinel(::Type{Tuple{Float64, Float64}}) = (0.0, 0.0)
@inline _run_key_sentinel(::Type{P}) where {P <: UnitSphericalPoint} = P(0.0, 0.0, 0.0)

#=
The run key as a struct rather than a 4-tuple, so its hash can stop at the node
pair for a run with no interior vertex — which is what a run between two
crossings of the same parent segment is, and on a crossing-dense curve that is
nearly all of them. A tuple key hashes all four fields unconditionally, and the
coordinate is two more mix rounds on a key that a single segment already
determines.
=#
struct _RunKey{P}
    lo   :: Int32
    hi   :: Int32
    nint :: Int32
    v    :: P
end

@inline Base.:(==)(a::_RunKey, b::_RunKey) =
    a.lo == b.lo && a.hi == b.hi && a.nint == b.nint && (a.nint == 0 || a.v == b.v)
@inline Base.isequal(a::_RunKey, b::_RunKey) = a == b

@inline function Base.hash(k::_RunKey, h::UInt)
    x = _nk_mix(UInt64(h), (UInt64(k.lo % UInt32) << 32) | UInt64(k.hi % UInt32))
    x = _nk_mix(x, UInt64(k.nint % UInt32))
    k.nint == 0 && return x % UInt
    return _nk_mix_point(x, k.v) % UInt
end

function _merge_noded_runs(arr::NodedArrangement{P}, sources::Vector{EdgeSourceInfo}) where {P}
    edgemap = Dict{_RunKey{P}, Int}()
    sizehint!(edgemap, length(arr.edges))
    merged = MergeEdge[]
    sizehint!(merged, length(arr.edges))
    sentinel = _run_key_sentinel(P)
    for ne in arr.edges
        src = sources[ne.string_idx]
        nint = ne.seg_hi - ne.seg_idx
        fwd = ne.node_lo <= ne.node_hi
        lo = fwd ? ne.node_lo : ne.node_hi
        hi = fwd ? ne.node_hi : ne.node_lo
        #-- the parent's coordinates are only read for a run that has an interior
        v = sentinel
        if nint != 0
            pts = arr.segstrings[ne.string_idx].pts
            v = fwd ? pts[ne.seg_idx + 1] : pts[ne.seg_hi]
        end
        key = _RunKey{P}(lo, hi, nint, v)
        idx = get(edgemap, key, 0)
        if idx == 0
            push!(merged, _merge_edge(ne, src))
            edgemap[key] = length(merged)
        else
            merged[idx] = _merge(merged[idx], _merge_edge(ne, src))
        end
    end
    return merged
end
