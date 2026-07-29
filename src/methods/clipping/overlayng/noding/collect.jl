# # Stage 1 — collect (design §2.3)
#
# Enumerate candidate segment pairs (A × B) through the reused RelateNG edge
# index and classify each with the exact kernel (`rk_classify_intersection`),
# recording symbolic node keys onto the parent segments. No intersection
# coordinate is ever constructed here.
#
# `seg_nodes` accumulates, per parent segment `(string_idx, seg_idx)`, the ids of
# nodes lying strictly in that segment's interior. `string_idx` is global into
# the arrangement's `segstrings`: A strings occupy `1:na`, B strings `na+1:end`.

@inline function _record_interior!(seg_nodes, string_idx::Int32, seg_idx::Int32, nid::Int32)
    push!(get!(() -> Int32[], seg_nodes, (string_idx, seg_idx)), nid)
    return nothing
end

function _collect_crossings!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ssa::AbstractVector{RelateSegmentString{P}},
        ssb::AbstractVector{RelateSegmentString{P}}, na::Int32;
        exact = True(), tree_a = nothing, tree_b = nothing) where {P}
    ta = tree_a === nothing ? _relate_edge_index(m, ssa) : tree_a
    tb = tree_b === nothing ? _relate_edge_index(m, ssb) : tree_b
    if !(ta === nothing || tb === nothing)
        SpatialTreeInterface.dual_depth_first_search(Extents.intersects, ta, tb) do ia, ib
            (sa, ka) = ta.data[ia]
            (sb, kb) = tb.data[ib]
            _classify_pair!(m, table, seg_nodes, ssa, sa, Int32(sa), Int32(ka),
                            ssb, sb, na + Int32(sb), Int32(kb); exact)
            return nothing
        end
    end
    #-- design §2.2 amendment: LINEAR inputs must additionally be self-noded
    _collect_self_crossings!(m, table, seg_nodes, ssa, Int32(0); exact)
    _collect_self_crossings!(m, table, seg_nodes, ssb, na; exact)
    return nothing
end

#=
Self-noding of one side's LINEAR segment strings.

The arrangement invariant is that no node lies strictly inside a noded edge
(design §2.1). The A×B pass alone does not establish it for linear inputs: a
node created by an A×B crossing can land in the interior of a *third* segment
belonging to the same input as one of the pair, and nothing then splits that
segment. `TestOverlayLA.xml` case 2 is the canonical instance — one
MULTILINESTRING component runs along a polygon hole boundary while a second
component of the same MULTILINESTRING crosses it, so the hole edge is split at
the two crossings but the collinear line component is not, the two no longer
share a node pair, and the edge merger never pairs them.

Only LINEAR strings are self-noded, and that is a cost decision, not a
correctness one. A valid linear input has no self-noding guarantee at all — a
MultiLineString may cross and even retrace itself — while a valid areal input
has almost one: its rings never cross or overlap, so the only self-incidence it
can carry is a single-point touch (a hole meeting its shell, or two
multipolygon components meeting). That touch point is always a vertex of one
ring, but it may lie in the *interior* of the other ring's segment, and then the
same gap opens for areas. That case is real and reproduced: it is the fuzz
suite's pinned `hole apex on the shell edge (self-touching input)` defect class,
which self-noding areal strings as well turns from 4 broken to 4 passing (the
one-line change is dropping the `DIM_L` filter below). It is not enabled because
it roughly doubles arrangement build time on real data (0.43 s -> 0.89 s over 28
Natural Earth 10m country pairs, 16k-68k points each) — the §2.2 trade — and
because a targeted version (each input's vertices against its own segments, the
only incidence a valid area can have) would very likely recover it for much
less. Polygon-only inputs skip the pass entirely: no index is built.
=#
function _collect_self_crossings!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ss::AbstractVector{RelateSegmentString{P}}, off::Int32; exact) where {P}
    lin = Int32[Int32(i) for i in eachindex(ss) if ss[i].dim == DIM_L]
    isempty(lin) && return nothing
    sub = [ss[i] for i in lin]
    #-- fewer than two segments in total: no pair to classify
    sum(s -> length(s.pts) - 1, sub; init = 0) < 2 && return nothing
    t = _relate_edge_index(m, sub)
    t === nothing && return nothing
    SpatialTreeInterface.dual_depth_first_search(Extents.intersects, t, t) do i1, i2
        (s1, k1) = t.data[i1]
        (s2, k2) = t.data[i2]
        #-- the dual traversal of a tree against itself yields every unordered
        #-- pair twice plus every self-pair; keep one ordering, drop self-pairs
        (s1 < s2 || (s1 == s2 && k1 < k2)) || return nothing
        _classify_pair!(m, table, seg_nodes, sub, s1, off + lin[s1], Int32(k1),
                        sub, s2, off + lin[s2], Int32(k2); exact)
        return nothing
    end
    return nothing
end

# Function barrier: statically-typed classification of one candidate pair
# (the do-blocks above are dynamic closures over the tree traversal). `gsa`/`gsb`
# are the *global* string indices (into the arrangement's `segstrings`) of the
# two segments' parents; `sa`/`sb` index the lists actually passed in.
function _classify_pair!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ssa, sa::Int, gsa::Int32, ksa::Int32,
        ssb, sb::Int, gsb::Int32, ksb::Int32; exact) where {P}
    a0 = ssa[sa].pts[ksa]; a1 = ssa[sa].pts[ksa + 1]
    b0 = ssb[sb].pts[ksb]; b1 = ssb[sb].pts[ksb + 1]

    cls = rk_classify_intersection(m, a0, a1, b0, b1; exact)
    kind = cls.kind
    if kind == SS_DISJOINT
        return nothing
    elseif kind == SS_PROPER
        #-- a proper crossing is strictly interior to BOTH segments
        nid = _intern_node!(table, crossing_node(a0, a1, b0, b1))
        _record_interior!(seg_nodes, gsa, ksa, nid)
        _record_interior!(seg_nodes, gsb, ksb, nid)
        return nothing
    end
    #-- SS_TOUCH / SS_COLLINEAR: every intersection point is an input vertex,
    #-- reported via the incidence flags (design §2.3, S1 census). The claim is
    #-- load-bearing; assert it.
    @assert (cls.a0_on_b || cls.a1_on_b || cls.b0_on_a || cls.b1_on_a) """
        touch/collinear intersection with no vertex incidence flag — \
        the touch=input-vertex invariant (design §2.3) is violated"""
    #-- a?/b? are vertices of A/B; record on the OTHER segment only where the
    #-- vertex lies strictly in that segment's interior (not at its endpoints).
    if cls.a0_on_b && a0 != b0 && a0 != b1
        _record_interior!(seg_nodes, gsb, ksb, _intern_node!(table, vertex_node(a0)))
    end
    if cls.a1_on_b && a1 != b0 && a1 != b1
        _record_interior!(seg_nodes, gsb, ksb, _intern_node!(table, vertex_node(a1)))
    end
    if cls.b0_on_a && b0 != a0 && b0 != a1
        _record_interior!(seg_nodes, gsa, ksa, _intern_node!(table, vertex_node(b0)))
    end
    if cls.b1_on_a && b1 != a0 && b1 != a1
        _record_interior!(seg_nodes, gsa, ksa, _intern_node!(table, vertex_node(b1)))
    end
    return nothing
end
