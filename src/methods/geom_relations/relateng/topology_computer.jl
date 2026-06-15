# # RelateNG topology computer
#
# Port of JTS `TopologyComputer.java` — the heart of the topology layer.
# Receives topological events from the evaluation phases (point locations,
# line ends, area vertices, edge intersections), translates each into
# DE-9IM dimension updates on the attached `TopologyPredicate`, and groups
# the `NodeSection`s of edge intersections per node for the final node
# topology analysis (`evaluate_nodes!`).
#
# Method order parallels the Java file, so this file diffs against its Java
# counterpart. Idiom changes:
#
# - The Java `Map<Coordinate, NodeSections> nodeMap` becomes a
#   `Dict{NodeKey, NodeSections}` keyed by the *symbolic* node identity
#   (design D2): vertex nodes key by their exact coordinate, proper-crossing
#   nodes by their canonicalized defining segment pair. No intersection
#   coordinate is ever constructed for node identity.
# - Because identical segment pairs always produce identical keys, the
#   "canonical example" of JTS's self-noding note (a self-crossing line
#   tested against a copy of one of the crossed segments) merges
#   automatically. Crossings of *different* segment pairs at the same
#   geometric point — which JTS only merges when floating-point intersection
#   coordinates happen to round identically — are merged *exactly* by the
#   D3 coincidence pass in `evaluate_nodes!` (self-noding predicates only).
# - The manifold and `exact` flag the kernel calls need are taken from
#   `geom_a` (the constructor asserts both inputs agree).

"""
    TopologyComputer(predicate, geom_a::RelateGeometry, geom_b::RelateGeometry)

The DE-9IM accumulation engine of RelateNG: translates topological events
into dimension updates on `predicate` and collects edge-intersection
[`NodeSection`](@ref)s per node (keyed by symbolic [`NodeKey`](@ref),
design D2) for `evaluate_nodes!`.

Port of JTS `TopologyComputer`.
"""
struct TopologyComputer{TP <: TopologyPredicate, RA <: RelateGeometry, RB <: RelateGeometry, P}
    predicate::TP
    geom_a::RA
    geom_b::RB
    node_sections::Dict{NodeKey{P}, NodeSections{P}}
end

function TopologyComputer(predicate::TopologyPredicate,
        geom_a::RelateGeometry, geom_b::RelateGeometry)
    #-- the kernel manifold/exact settings are read from geom_a below;
    #-- both inputs must have been built with the same settings
    (geom_a.m == geom_b.m && geom_a.exact == geom_b.exact) ||
        throw(ArgumentError("RelateGeometry manifold/exact settings of the A and B inputs must agree"))
    #-- P is the manifold's kernel point type, matching the coordinate type of
    #-- every segment string / node point produced at ingest (Tuple{Float64,
    #-- Float64} for Planar, UnitSphericalPoint{Float64} for Spherical)
    P = _kernel_point_type(geom_a.m)
    tc = TopologyComputer(predicate, geom_a, geom_b, Dict{NodeKey{P}, NodeSections{P}}())
    init_exterior_dims!(tc)
    return tc
end

# The manifold / exactness flag for kernel calls (asserted equal across both
# inputs in the constructor).
_manifold(tc::TopologyComputer) = tc.geom_a.m
_exact(tc::TopologyComputer) = tc.geom_a.exact

# Port of TopologyComputer.initExteriorDims (private): determine a priori
# partial EXTERIOR topology based on the real dimensions.
function init_exterior_dims!(tc::TopologyComputer)
    dim_real_a = get_dimension_real(tc.geom_a)
    dim_real_b = get_dimension_real(tc.geom_b)

    if dim_real_a == DIM_P && dim_real_b == DIM_L
        #-- For P/L case, P exterior intersects L interior
        update_dim!(tc, LOC_EXTERIOR, LOC_INTERIOR, DIM_L)
    elseif dim_real_a == DIM_L && dim_real_b == DIM_P
        update_dim!(tc, LOC_INTERIOR, LOC_EXTERIOR, DIM_L)
    elseif dim_real_a == DIM_P && dim_real_b == DIM_A
        #-- For P/A case, the Area Int and Bdy intersect the Point exterior.
        update_dim!(tc, LOC_EXTERIOR, LOC_INTERIOR, DIM_A)
        update_dim!(tc, LOC_EXTERIOR, LOC_BOUNDARY, DIM_L)
    elseif dim_real_a == DIM_A && dim_real_b == DIM_P
        update_dim!(tc, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
        update_dim!(tc, LOC_BOUNDARY, LOC_EXTERIOR, DIM_L)
    elseif dim_real_a == DIM_L && dim_real_b == DIM_A
        update_dim!(tc, LOC_EXTERIOR, LOC_INTERIOR, DIM_A)
    elseif dim_real_a == DIM_A && dim_real_b == DIM_L
        update_dim!(tc, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
    elseif dim_real_a == DIM_FALSE || dim_real_b == DIM_FALSE
        #-- cases where one geom is EMPTY
        if dim_real_a != DIM_FALSE
            init_exterior_empty!(tc, GEOM_A)
        end
        if dim_real_b != DIM_FALSE
            init_exterior_empty!(tc, GEOM_B)
        end
    end
    return nothing
end

# Port of TopologyComputer.initExteriorEmpty (private).
function init_exterior_empty!(tc::TopologyComputer, geom_non_empty::Bool)
    dim_non_empty = get_dimension(tc, geom_non_empty)
    if dim_non_empty == DIM_P
        update_dim!(tc, geom_non_empty, LOC_INTERIOR, LOC_EXTERIOR, DIM_P)
    elseif dim_non_empty == DIM_L
        if has_boundary(get_geometry(tc, geom_non_empty))
            update_dim!(tc, geom_non_empty, LOC_BOUNDARY, LOC_EXTERIOR, DIM_P)
        end
        update_dim!(tc, geom_non_empty, LOC_INTERIOR, LOC_EXTERIOR, DIM_L)
    elseif dim_non_empty == DIM_A
        update_dim!(tc, geom_non_empty, LOC_BOUNDARY, LOC_EXTERIOR, DIM_L)
        update_dim!(tc, geom_non_empty, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
    end
    return nothing
end

# Port of TopologyComputer.getGeometry (private).
get_geometry(tc::TopologyComputer, is_a::Bool) = is_a ? tc.geom_a : tc.geom_b

# Port of TopologyComputer.getDimension.
get_dimension(tc::TopologyComputer, is_a::Bool) = get_dimension(get_geometry(tc, is_a))

# Port of TopologyComputer.isAreaArea.
is_area_area(tc::TopologyComputer) =
    get_dimension(tc, GEOM_A) == DIM_A && get_dimension(tc, GEOM_B) == DIM_A

"""
    is_self_noding_required(tc::TopologyComputer)

Indicates whether the input geometries require self-noding for correct
evaluation of specific spatial predicates. Self-noding is required for
geometries which may have self-crossing linework, or may have lines lying
in the boundary of an area. This ensures that node locations match in
situations where a self-crossing and mutual crossing occur at the same
logical location (here via the D3 coincidence-merge pass, since node
identities are symbolic).

Currently self-noding is required for:
- A geoms which require self-noding (lines or GCs, except for single-polygon GCs)
- B geoms which are mixed A/L GCs

Port of TopologyComputer.isSelfNodingRequired.
"""
function is_self_noding_required(tc::TopologyComputer)
    require_self_noding(tc.predicate) || return false

    is_self_noding_required(tc.geom_a) && return true

    #-- if B is a mixed GC with A and L require full noding
    has_area_and_line(tc.geom_b) && return true

    return false
end

# Port of TopologyComputer.isExteriorCheckRequired.
is_exterior_check_required(tc::TopologyComputer, is_a::Bool) =
    require_exterior_check(tc.predicate, is_a)

# Port of TopologyComputer.updateDim(locA, locB, dimension) (private).
function update_dim!(tc::TopologyComputer, loc_a::Integer, loc_b::Integer, dim::Integer)
    update_dim!(tc.predicate, loc_a, loc_b, dim)
    return nothing
end

# Port of TopologyComputer.updateDim(isAB, loc1, loc2, dimension) (private):
# `loc1`/`loc2` are ordered source/target; swapped when the source is B.
function update_dim!(tc::TopologyComputer, is_ab::Bool, loc1::Integer, loc2::Integer, dim::Integer)
    if is_ab
        update_dim!(tc, loc1, loc2, dim)
    else
        #-- is ordered BA
        update_dim!(tc, loc2, loc1, dim)
    end
    return nothing
end

# Port of TopologyComputer.isResultKnown.
is_result_known(tc::TopologyComputer) = is_known(tc.predicate)

# Port of TopologyComputer.getResult.
get_result(tc::TopologyComputer) = predicate_value(tc.predicate)

# Port of TopologyComputer.finish: finalize the evaluation.
function finish!(tc::TopologyComputer)
    finish!(tc.predicate)
    return nothing
end

# Port of TopologyComputer.getNodeSections (private); the map is keyed by
# the symbolic NodeKey instead of a Coordinate (design D2).
_get_node_sections(tc::TopologyComputer, node::NodeKey) =
    get!(() -> NodeSections(node), tc.node_sections, node)

# Port of TopologyComputer.addIntersection.
function add_intersection!(tc::TopologyComputer, a::NodeSection, b::NodeSection)
    if !is_same_geometry(a, b)
        update_intersection_ab!(tc, a, b)
    end
    #-- add edges to node to allow full topology evaluation later
    add_node_sections!(tc, a, b)
    return nothing
end

# Port of TopologyComputer.updateIntersectionAB (private): update topology
# for an intersection between A and B.
function update_intersection_ab!(tc::TopologyComputer, a::NodeSection, b::NodeSection)
    if is_area_area(a, b)
        update_area_area_cross!(tc, a, b)
    end
    update_node_location!(tc, a, b)
    return nothing
end

#=
Port of TopologyComputer.updateAreaAreaCross (private): updates topology for
an AB Area-Area crossing node. Sections cross at a node if (a) the
intersection is proper (i.e. in the interior of two segments) or (b) if
non-proper then whether the linework crosses is determined by the geometry
of the segments on either side of the node. In these situations the area
geometry interiors intersect (in dimension 2).

A proper intersection short-circuits `rk_is_crossing` (which requires a
vertex-node apex; proper crossings cross by construction).
=#
function update_area_area_cross!(tc::TopologyComputer, a::NodeSection, b::NodeSection)
    if is_proper(a, b) || rk_is_crossing(_manifold(tc), node_pt(a),
            get_vertex(a, 0), get_vertex(a, 1),
            get_vertex(b, 0), get_vertex(b, 1); exact = _exact(tc))
        update_dim!(tc, LOC_INTERIOR, LOC_INTERIOR, DIM_A)
    end
    return nothing
end

# Port of TopologyComputer.updateNodeLocation (private): updates topology
# for a node at an AB edge intersection. The Java passes the node
# Coordinate; here the symbolic NodeKey is located (see `locate_node` on
# NodeKey below).
function update_node_location!(tc::TopologyComputer, a::NodeSection, b::NodeSection)
    pt = node_pt(a)
    loc_a = locate_node(tc.geom_a, pt, get_polygonal(a))
    loc_b = locate_node(tc.geom_b, pt, get_polygonal(b))
    update_dim!(tc, loc_a, loc_b, DIM_P)
    return nothing
end

# Port of TopologyComputer.addNodeSections (private).
function add_node_sections!(tc::TopologyComputer, ns0::NodeSection, ns1::NodeSection)
    sections = _get_node_sections(tc, node_pt(ns0))
    add_node_section!(sections, ns0)
    add_node_section!(sections, ns1)
    return nothing
end

# Port of TopologyComputer.addPointOnPointInterior.
function add_point_on_point_interior!(tc::TopologyComputer, pt)
    update_dim!(tc, LOC_INTERIOR, LOC_INTERIOR, DIM_P)
    return nothing
end

# Port of TopologyComputer.addPointOnPointExterior.
function add_point_on_point_exterior!(tc::TopologyComputer, is_geom_a::Bool, pt)
    update_dim!(tc, is_geom_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_P)
    return nothing
end

# Port of TopologyComputer.addPointOnGeometry.
function add_point_on_geometry!(tc::TopologyComputer, is_point_a::Bool,
        loc_target::Integer, dim_target::Integer, pt)
    #-- update entry for Point interior
    update_dim!(tc, is_point_a, LOC_INTERIOR, loc_target, DIM_P)

    #-- an empty geometry has no points to infer entries from
    is_geom_empty(get_geometry(tc, !is_point_a)) && return nothing

    if dim_target == DIM_P
        return nothing
    elseif dim_target == DIM_L
        #=
        Because zero-length lines are handled, a point lying in the exterior
        of the line target may imply either P or L for the Exterior
        interaction
        =#
        #TODO: determine if effective dimension of linear target is L?
        return nothing
    elseif dim_target == DIM_A
        #=
        If a point intersects an area target, then the area interior and
        boundary must extend beyond the point and thus interact with its
        exterior.
        =#
        update_dim!(tc, is_point_a, LOC_EXTERIOR, LOC_INTERIOR, DIM_A)
        update_dim!(tc, is_point_a, LOC_EXTERIOR, LOC_BOUNDARY, DIM_L)
        return nothing
    end
    throw(ArgumentError("unknown target dimension: $dim_target"))
end

"""
    add_line_end_on_geometry!(tc, is_line_a, loc_line_end, loc_target, dim_target, pt)

Add topology for a line end. The line end point must be "significant"; i.e.
not contained in an area if the source is a mixed-dimension GC.
`loc_line_end` is the location of the line end (Interior or Boundary);
`loc_target` the location on the target geometry; `dim_target` the dimension
of the interacting target geometry element (if any), or the dimension of
the target.

Port of TopologyComputer.addLineEndOnGeometry.
"""
function add_line_end_on_geometry!(tc::TopologyComputer, is_line_a::Bool,
        loc_line_end::Integer, loc_target::Integer, dim_target::Integer, pt)
    #-- record topology at line end point
    update_dim!(tc, is_line_a, loc_line_end, loc_target, DIM_P)

    #-- an empty geometry has no points to infer entries from
    is_geom_empty(get_geometry(tc, !is_line_a)) && return nothing

    #-- Line and Area targets may have additional topology
    if dim_target == DIM_P
        return nothing
    elseif dim_target == DIM_L
        add_line_end_on_line!(tc, is_line_a, loc_line_end, loc_target, pt)
        return nothing
    elseif dim_target == DIM_A
        add_line_end_on_area!(tc, is_line_a, loc_line_end, loc_target, pt)
        return nothing
    end
    throw(ArgumentError("unknown target dimension: $dim_target"))
end

# Port of TopologyComputer.addLineEndOnLine (private).
function add_line_end_on_line!(tc::TopologyComputer, is_line_a::Bool,
        loc_line_end::Integer, loc_line::Integer, pt)
    #=
    When a line end is in the EXTERIOR of a Line, some length of the source
    Line INTERIOR is also in the target Line EXTERIOR. This works for
    zero-length lines as well.
    =#
    if loc_line == LOC_EXTERIOR
        update_dim!(tc, is_line_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_L)
    end
    return nothing
end

# Port of TopologyComputer.addLineEndOnArea (private).
function add_line_end_on_area!(tc::TopologyComputer, is_line_a::Bool,
        loc_line_end::Integer, loc_area::Integer, pt)
    if loc_area != LOC_BOUNDARY
        #=
        When a line end is in an Area INTERIOR or EXTERIOR some length of
        the source Line Interior AND the Exterior of the line is also in
        that location of the target.
        NOTE: this assumes the line end is NOT also in an Area of a mixed-dim GC
        =#
        #TODO: handle zero-length lines?
        update_dim!(tc, is_line_a, LOC_INTERIOR, loc_area, DIM_L)
        update_dim!(tc, is_line_a, LOC_EXTERIOR, loc_area, DIM_A)
    end
    return nothing
end

"""
    add_area_vertex!(tc, is_area_a, loc_area, loc_target, dim_target, pt)

Adds topology for an area vertex interaction with a target geometry element.
Assumes the target geometry element has highest dimension (i.e. if the point
lies on two elements of different dimension, the location on the higher
dimension element is provided. This is the semantic provided by
[`RelatePointLocator`](@ref).)

Note that in a GeometryCollection containing overlapping or adjacent
polygons, the area vertex location may be INTERIOR instead of BOUNDARY.

Port of TopologyComputer.addAreaVertex.
"""
function add_area_vertex!(tc::TopologyComputer, is_area_a::Bool,
        loc_area::Integer, loc_target::Integer, dim_target::Integer, pt)
    if loc_target == LOC_EXTERIOR
        update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
        #=
        If area vertex is on Boundary further topology can be deduced from
        the neighbourhood around the boundary vertex. This is always the
        case for polygonal geometries. For GCs, the vertex may be either on
        boundary or in interior (i.e. of overlapping or adjacent polygons)
        =#
        if loc_area == LOC_BOUNDARY
            update_dim!(tc, is_area_a, LOC_BOUNDARY, LOC_EXTERIOR, DIM_L)
            update_dim!(tc, is_area_a, LOC_EXTERIOR, LOC_EXTERIOR, DIM_A)
        end
        return nothing
    end
    if dim_target == DIM_P
        add_area_vertex_on_point!(tc, is_area_a, loc_area, pt)
        return nothing
    elseif dim_target == DIM_L
        add_area_vertex_on_line!(tc, is_area_a, loc_area, loc_target, pt)
        return nothing
    elseif dim_target == DIM_A
        add_area_vertex_on_area!(tc, is_area_a, loc_area, loc_target, pt)
        return nothing
    end
    throw(ArgumentError("unknown target dimension: $dim_target"))
end

#=
Port of TopologyComputer.addAreaVertexOnPoint (private): updates topology
for an area vertex (in Interior or on Boundary) intersecting a point. Note
that because the largest dimension of intersecting target is determined,
the intersecting point is not part of any other target geometry, and hence
its neighbourhood is in the Exterior of the target.
=#
function add_area_vertex_on_point!(tc::TopologyComputer, is_area_a::Bool,
        loc_area::Integer, pt)
    #-- Assert: loc_area != EXTERIOR
    #-- Assert: loc_target == INTERIOR
    #-- The vertex location intersects the Point.
    update_dim!(tc, is_area_a, loc_area, LOC_INTERIOR, DIM_P)
    #-- The area interior intersects the point's exterior neighbourhood.
    update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
    #=
    If the area vertex is on the boundary, the area boundary and exterior
    intersect the point's exterior neighbourhood
    =#
    if loc_area == LOC_BOUNDARY
        update_dim!(tc, is_area_a, LOC_BOUNDARY, LOC_EXTERIOR, DIM_L)
        update_dim!(tc, is_area_a, LOC_EXTERIOR, LOC_EXTERIOR, DIM_A)
    end
    return nothing
end

# Port of TopologyComputer.addAreaVertexOnLine (private).
function add_area_vertex_on_line!(tc::TopologyComputer, is_area_a::Bool,
        loc_area::Integer, loc_target::Integer, pt)
    #-- Assert: loc_area != EXTERIOR
    #=
    If an area vertex intersects a line, all we know is the intersection at
    that point. e.g. the line may or may not be collinear with the area
    boundary, and the line may or may not intersect the area interior.
    Full topology is determined later by node analysis
    =#
    update_dim!(tc, is_area_a, loc_area, loc_target, DIM_P)
    if loc_area == LOC_INTERIOR
        #-- The area interior intersects the line's exterior neighbourhood.
        update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
    end
    return nothing
end

# Port of TopologyComputer.addAreaVertexOnArea (public in Java).
function add_area_vertex_on_area!(tc::TopologyComputer, is_area_a::Bool,
        loc_area::Integer, loc_target::Integer, pt)
    if loc_target == LOC_BOUNDARY
        if loc_area == LOC_BOUNDARY
            #-- B/B topology is fully computed later by node analysis
            update_dim!(tc, is_area_a, LOC_BOUNDARY, LOC_BOUNDARY, DIM_P)
        else
            #-- loc_area == INTERIOR
            update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_INTERIOR, DIM_A)
            update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_BOUNDARY, DIM_L)
            update_dim!(tc, is_area_a, LOC_INTERIOR, LOC_EXTERIOR, DIM_A)
        end
    else
        #-- loc_target is INTERIOR or EXTERIOR
        update_dim!(tc, is_area_a, LOC_INTERIOR, loc_target, DIM_A)
        #=
        If area vertex is on Boundary further topology can be deduced from
        the neighbourhood around the boundary vertex. This is always the
        case for polygonal geometries. For GCs, the vertex may be either on
        boundary or in interior (i.e. of overlapping or adjacent polygons)
        =#
        if loc_area == LOC_BOUNDARY
            update_dim!(tc, is_area_a, LOC_BOUNDARY, loc_target, DIM_L)
            update_dim!(tc, is_area_a, LOC_EXTERIOR, loc_target, DIM_A)
        end
    end
    return nothing
end

# Port of TopologyComputer.evaluateNodes, preceded by the D3
# coincidence-merge pass (not present in Java, where concrete node
# coordinates merge in the HashMap when they round identically).
function evaluate_nodes!(tc::TopologyComputer)
    _merge_coincident_nodes!(tc)
    for node_sections in values(tc.node_sections)
        if has_interaction_ab(node_sections)
            evaluate_node!(tc, node_sections)
            is_result_known(tc) && return nothing
        end
    end
    return nothing
end

#=
D3 coincidence-merge pass: distinct symbolic crossing keys (different
segment pairs) — and vertex keys — may denote the same geometric point.
Group them exactly and merge their NodeSections into one node.

In JTS this merging is implicit and unconditional: the nodeMap is keyed by
the *constructed* intersection Coordinate, so a proper crossing whose
(floating-point) intersection point coincides with a vertex node — or with
another crossing — lands in the same map entry, in every mode (not only
under self-noding; e.g. RelateNGTest.testPolygonLineCrossingContained needs
a B-line proper crossing of one A polygon merged with its vertex touch of
another). Here node identity is symbolic, so the merge is an explicit pass,
run whenever any crossing key exists.

Candidate grouping is by *exact bounding boxes* (the F1 follow-up to
design D3): a vertex key's box is its exact coordinate; a crossing key's
box is the intersection of its two defining segments' bounding boxes,
which provably contains the exact crossing point (the point lies on both
segments). Boxes are computed with exact Float64 comparisons — no
rounding is involved — so coincident keys (equal exact points) always
have overlapping boxes and always land in the same overlap cluster.
Clusters are the sweep closure of box overlap along x, then y: a
conservative superset of the true coincidence classes. Within a
multi-member cluster containing at least one crossing key, coincidence is
confirmed *exactly* via the rational `_exact_node_point` (distinct exact
points whose boxes happen to overlap are never merged); vertex-only
clusters need no check at all, since distinct vertex keys key by their
exact coordinates and can never coincide. This keeps the merge decisions
exact (design D3) at sorting cost O(N log N), with the rational
arithmetic reserved for genuinely near-coincident nodes instead of every
crossing key on every evaluation.

A vertex key is preferred as the canonical merged node: its coordinate is
exact, so the edge wheel and node location never need the rational apex.
Otherwise the merged crossing node's wheel compares foreign directions
around the exact rational apex (`rk_compare_edge_dir` slow path).
=#
function _merge_coincident_nodes!(tc::TopologyComputer)
    nodemap = tc.node_sections
    length(nodemap) > 1 || return nothing
    any(k -> k.is_crossing, keys(nodemap)) || return nothing
    K = keytype(nodemap)
    #-- collect (x interval, y interval, key) and cluster by box overlap
    items = Vector{Tuple{NTuple{2, Float64}, NTuple{2, Float64}, K}}()
    sizehint!(items, length(nodemap))
    for k in keys(nodemap)
        xint, yint = _node_key_box(k)
        push!(items, (xint, yint, k))
    end
    sort!(items; by = it -> it[1][1])
    i = 1
    n = length(items)
    while i <= n
        #-- extend the x-cluster while the next box starts inside it
        j = i
        xhi = items[i][1][2]
        while j < n && items[j + 1][1][1] <= xhi
            j += 1
            xhi = max(xhi, items[j][1][2])
        end
        j > i && _merge_coincident_y_clusters!(nodemap, items[i:j])
        i = j + 1
    end
    return nothing
end

# The exact coordinate box guaranteed to contain a node key's point: the
# vertex coordinate itself, or the intersection of the defining segments'
# bounding boxes for a proper crossing.
function _node_key_box(k::NodeKey)
    if !k.is_crossing
        return (k.pt[1], k.pt[1]), (k.pt[2], k.pt[2])
    end
    axlo, axhi = minmax(k.pt[1], k.a1[1])
    aylo, ayhi = minmax(k.pt[2], k.a1[2])
    bxlo, bxhi = minmax(k.b0[1], k.b1[1])
    bylo, byhi = minmax(k.b0[2], k.b1[2])
    return (max(axlo, bxlo), min(axhi, bxhi)), (max(aylo, bylo), min(ayhi, byhi))
end

# Second sweep dimension: within one x-overlap cluster, cluster by y
# overlap; multi-member y-clusters with a crossing key are the candidate
# coincidence groups handed to exact confirmation.
function _merge_coincident_y_clusters!(nodemap::Dict, items::Vector)
    sort!(items; by = it -> it[2][1])
    i = 1
    n = length(items)
    while i <= n
        j = i
        yhi = items[i][2][2]
        while j < n && items[j + 1][2][1] <= yhi
            j += 1
            yhi = max(yhi, items[j][2][2])
        end
        if j > i
            group = [items[t][3] for t in i:j]
            any(k -> k.is_crossing, group) &&
                _merge_coincident_group!(nodemap, group)
        end
        i = j + 1
    end
    return nothing
end

# Merge each exact-coincidence class within one rounded-coordinate bucket.
function _merge_coincident_group!(nodemap::Dict, group::Vector)
    #-- confirm coincidence exactly (rational arithmetic): bucket members
    #-- are only candidates, since distinct exact points may round together
    exacts = [_exact_node_point(k) for k in group]
    merged = falses(length(group))
    for i in eachindex(group)
        merged[i] && continue
        cls = [i]
        for j in (i + 1):length(group)
            merged[j] && continue
            if exacts[j] == exacts[i]
                push!(cls, j)
                merged[j] = true
            end
        end
        length(cls) > 1 || continue
        #-- prefer a vertex key as the canonical node (at most one exists:
        #-- distinct vertex keys are distinct points)
        ci = findfirst(j -> !group[j].is_crossing, cls)
        canonical = group[cls[ci === nothing ? 1 : ci]]
        target = nodemap[canonical]
        for j in cls
            k = group[j]
            k == canonical && continue
            _merge_node_sections!(target, pop!(nodemap, k), canonical)
        end
    end
    return nothing
end

# Move the sections of `src` into `target`, rewriting their node to the
# canonical key so the section angle comparators and `create_node` see a
# consistent apex.
function _merge_node_sections!(target::NodeSections, src::NodeSections, canonical::NodeKey)
    for ns in src.sections
        add_node_section!(target, _section_with_node(ns, canonical))
    end
    return nothing
end

_section_with_node(ns::NodeSection, node::NodeKey) =
    NodeSection(ns.is_a, ns.dim, ns.id, ns.ring_id, ns.polygonal,
        ns.is_node_at_vertex, ns.v0, node, ns.v1)

# Port of TopologyComputer.evaluateNode (private).
function evaluate_node!(tc::TopologyComputer, node_sections::NodeSections)
    p = get_coordinate(node_sections)
    node = create_node(_manifold(tc), node_sections; exact = _exact(tc))
    #-- Node must have edges for geom, but may also be in interior of a overlapping GC
    is_area_interior_a = is_node_in_area(tc.geom_a, p, get_polygonal(node_sections, GEOM_A))
    is_area_interior_b = is_node_in_area(tc.geom_b, p, get_polygonal(node_sections, GEOM_B))
    finish!(node, is_area_interior_a, is_area_interior_b)
    evaluate_node_edges!(tc, node)
    return nothing
end

# Port of TopologyComputer.evaluateNodeEdges (private).
function evaluate_node_edges!(tc::TopologyComputer, node::RelateNode)
    #TODO: collect distinct dim settings by using temporary matrix?
    for e in get_edges(node)
        #-- An optimization to avoid updates for cases with a linear geometry
        if is_area_area(tc)
            update_dim!(tc, edge_location(e, GEOM_A, POS_LEFT),
                edge_location(e, GEOM_B, POS_LEFT), DIM_A)
            update_dim!(tc, edge_location(e, GEOM_A, POS_RIGHT),
                edge_location(e, GEOM_B, POS_RIGHT), DIM_A)
        end
        update_dim!(tc, edge_location(e, GEOM_A, POS_ON),
            edge_location(e, GEOM_B, POS_ON), DIM_L)
    end
    return nothing
end

#==========================================================================
## Locating symbolic nodes in an input geometry

Java's updateNodeLocation/evaluateNode pass the node Coordinate into
RelateGeometry.locateNode/isNodeInArea. Here the node may be a symbolic
proper-crossing key (design D2) with no stored coordinate:

- Vertex keys locate by their exact coordinate, as in Java.
- Crossing keys on a *polygonal* geometry need no coordinate at all: a
  node of a Polygon/MultiPolygon lies on its boundary (the same exact
  shortcut as `locate_with_dim`'s isNode branch).
- Otherwise (lineal geometries and GCs, where another element may cover
  the node) a representative coordinate is required. The exact rational
  crossing point is computed and rounded to Float64 — at least as precise
  as JTS, whose node coordinate is the floating-point intersection
  computed by RobustLineIntersector.
==========================================================================#

function locate_node(rg::RelateGeometry, key::NodeKey, parent_polygonal)
    key.is_crossing || return locate_node(rg, key.pt, parent_polygonal)
    #-- exact shortcut: a node of a polygonal geometry is on its boundary
    if GI.trait(rg.geom) isa Union{GI.PolygonTrait, GI.MultiPolygonTrait}
        return LOC_BOUNDARY
    end
    return locate_node(rg, _crossing_locate_point(key), parent_polygonal)
end

function is_node_in_area(rg::RelateGeometry, key::NodeKey, parent_polygonal)
    key.is_crossing || return is_node_in_area(rg, key.pt, parent_polygonal)
    #-- exact shortcut: a node of a polygonal geometry is on its boundary,
    #-- never in its interior
    if GI.trait(rg.geom) isa Union{GI.PolygonTrait, GI.MultiPolygonTrait}
        return false
    end
    return is_node_in_area(rg, _crossing_locate_point(key), parent_polygonal)
end

# Representative Float64 coordinate of a proper-crossing node,
# deterministically rounded via 256-bit BigFloat from the exact rational
# crossing point (BigFloat avoids overflow in the Rational → Float64
# conversion of huge numerators/denominators; the double rounding through
# BigFloat's default precision is deterministic, though not strictly
# correctly rounded).
function _crossing_locate_point(key::NodeKey)
    xr, yr = _exact_crossing_point(key.pt, key.a1, key.b0, key.b1)
    return (Float64(BigFloat(xr)), Float64(BigFloat(yr)))
end
