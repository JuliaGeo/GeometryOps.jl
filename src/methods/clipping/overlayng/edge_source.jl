# # EdgeSourceInfo — per-ring/line source topology (port of JTS `EdgeSourceInfo`)
#
# Phase 2a of the OverlayNG port (design doc §2.7, §3 amendment 3). Records the
# topological information carried by one source linework (a polygon ring or a
# line) from a single input geometry: which input it came from, its dimension,
# its ring role, and — for area rings — the signed `depth_delta`.
#
# One `EdgeSourceInfo` exists per `RelateSegmentString`, indexed by the segment
# string's position in the arrangement's `segstrings`. A `NodedEdge` reaches its
# source info through `NodedEdge.string_idx`; nothing is copied per noded edge
# (design §2.1 invariant 5), so nothing can desynchronize.
#
# Everything here is internal to GeometryOps — nothing is exported.

"""
    EdgeSourceInfo

Source topology of one input linework (port of JTS `EdgeSourceInfo`): `index`
(`0` = A, `1` = B), `dim` (`DIM_A` for a ring, `DIM_L` for a line), `is_hole`
(ring role), and `depth_delta` (the signed side-labelling delta of an area ring,
`0` for a line). Consumed by the edge merger to build `OverlayLabel`s.
"""
struct EdgeSourceInfo
    index       :: Int8
    dim         :: Int8
    is_hole     :: Bool
    depth_delta :: Int8
end

# The depth delta of one area ring, derived ONCE from the material-interior
# authority (design §2.7). JTS `EdgeNodingBuilder.computeDepthDelta` assigns the
# canonical delta `+1` to a ring in "Exterior-on-Left" orientation and flips to
# `-1` otherwise; equivalently (`Edge.locationLeft/Right`) a positive delta means
# Left = EXTERIOR, Right = INTERIOR — i.e. the material interior lies on the
# RIGHT. So `depth_delta = material_interior_on_left ? -1 : +1`, computed on the
# ring's stored order and matching that convention exactly. Folding the hole flip
# into `_ring_material_interior_on_left` keeps relate / overlay / extents in
# agreement by construction.
function _ring_depth_delta(m::Manifold, pts::Vector, is_hole::Bool; exact)
    return _ring_material_interior_on_left(m, pts, is_hole; exact) ? Int8(-1) : Int8(1)
end

# Build the source info for one segment string. Area rings carry a depth delta;
# lines carry none (`dim == DIM_L`, `depth_delta == 0`).
function _edge_source_info(m::Manifold, ss::RelateSegmentString; exact = True())
    index = ss.is_a ? Int8(0) : Int8(1)
    if ss.dim == DIM_A
        is_hole = _ss_is_hole(ss)
        return EdgeSourceInfo(index, DIM_A, is_hole, _ring_depth_delta(m, ss.pts, is_hole; exact))
    end
    #-- line (or any non-area edge): no side labelling
    return EdgeSourceInfo(index, DIM_L, false, Int8(0))
end

# The per-string source-info table of an arrangement, aligned to `segstrings`
# (so `sources[NodedEdge.string_idx]` is the noded edge's source info).
function _edge_source_infos(m::Manifold, arr::NodedArrangement; exact = True())
    return [_edge_source_info(m, ss; exact) for ss in arr.segstrings]
end
