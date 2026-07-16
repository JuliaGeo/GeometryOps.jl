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
    (ta === nothing || tb === nothing) && return nothing
    SpatialTreeInterface.dual_depth_first_search(Extents.intersects, ta, tb) do ia, ib
        (sa, ka) = ta.data[ia]
        (sb, kb) = tb.data[ib]
        _classify_pair!(m, table, seg_nodes, ssa, ssb, na, sa, ka, sb, kb; exact)
        return nothing
    end
    return nothing
end

# Function barrier: statically-typed classification of one candidate pair
# (the do-block above is a dynamic closure over the tree traversal).
function _classify_pair!(m::Manifold, table::NodeTable{P}, seg_nodes,
        ssa, ssb, na::Int32, sa::Int, ka::Int, sb::Int, kb::Int; exact) where {P}
    a0 = ssa[sa].pts[ka]; a1 = ssa[sa].pts[ka + 1]
    b0 = ssb[sb].pts[kb]; b1 = ssb[sb].pts[kb + 1]
    gsa = Int32(sa)                # A string global index
    gsb = na + Int32(sb)           # B string global index
    ksa = Int32(ka); ksb = Int32(kb)

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
