# Tests for the RelateNG node-topology layer (node_sections.jl and
# polygon_node_converter.jl): `NodeSection` accessors and comparators, the
# `NodeSections` collector, and the `PolygonNodeConverter` port.
# Ports of JTS NodeSection.java / NodeSections.java / PolygonNodeConverter.java.
# No dedicated JUnit file exists for the first two classes; those tests follow
# the implementation plan (Task 15): EdgeAngleComparator ordering,
# isProper/isNodeAtVertex semantics, the prepareSections ordering invariant,
# and the (partial, see node_sections.jl) createNode port on a simple
# two-area touch. The converter tests port PolygonNodeConverterTest.java
# in full (Task 16).

using Test
import GeometryOps as GO
import GeometryOps: Planar, True
import GeoInterface as GI

const NODE = GO.vertex_node((0.0, 0.0))

# Convenience constructor mirroring the Java NodeSection constructor
# argument order, with defaults for an A-geometry shell corner at NODE.
make_section(; is_a = true, dim = GO.DIM_A, id = 1, ring_id = 0, poly = nothing,
        at_vertex = true, v0 = (1.0, 0.0), node = NODE, v1 = (0.0, 1.0)) =
    GO.NodeSection(is_a, Int8(dim), Int32(id), Int32(ring_id), poly, at_vertex, v0, node, v1)

@testset "NodeSection accessors" begin
    poly = GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 0.0)]])
    ns = make_section(; poly)
    @test GO.get_vertex(ns, 0) == (1.0, 0.0)
    @test GO.get_vertex(ns, 1) == (0.0, 1.0)
    @test GO.node_pt(ns) === NODE
    @test GO.dimension(ns) == GO.DIM_A
    @test GO.id(ns) == 1
    @test GO.ring_id(ns) == 0
    @test GO.get_polygonal(ns) === poly
    @test GO.is_shell(ns)
    @test GO.is_area(ns)
    @test GO.is_a(ns)

    hole = make_section(; ring_id = 1)
    @test !GO.is_shell(hole)

    line = make_section(; is_a = false, dim = GO.DIM_L, id = 2, ring_id = -1,
        at_vertex = false, v1 = nothing)
    @test !GO.is_area(line)
    @test !GO.is_a(line)
    @test GO.get_polygonal(line) === nothing
    @test GO.get_vertex(line, 1) === nothing

    @test GO.is_area_area(ns, hole)
    @test !GO.is_area_area(ns, line)

    @test GO.is_same_geometry(ns, hole)
    @test !GO.is_same_geometry(ns, line)
    @test GO.is_same_polygon(ns, hole)              # same geometry, same id
    @test !GO.is_same_polygon(ns, make_section(; id = 2))
    @test !GO.is_same_polygon(ns, make_section(; is_a = false))

    # toString port (smoke test only)
    @test sprint(show, ns) isa String
    @test sprint(show, line) isa String
end

@testset "is_proper / is_node_at_vertex" begin
    at_v = make_section(; at_vertex = true)
    proper = make_section(; at_vertex = false)
    @test GO.is_node_at_vertex(at_v)
    @test !GO.is_node_at_vertex(proper)
    # Java isProper() == !isNodeAtVertex(): "proper" means the node falls in
    # a segment interior, not at a vertex of the section's component.
    @test !GO.is_proper(at_v)
    @test GO.is_proper(proper)
    # static isProper(a, b): both sections proper
    @test GO.is_proper(proper, proper)
    @test !GO.is_proper(proper, at_v)
    @test !GO.is_proper(at_v, at_v)
end

@testset "EdgeAngleComparator" begin
    m = Planar()
    # Directions in CCW angular order from the positive X-axis around (0,0).
    dirs_ccw = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0),
        (-1.0, 0.0), (-1.0, -1.0), (0.0, -1.0), (1.0, -1.0)]
    shuffled = dirs_ccw[[5, 1, 8, 3, 7, 2, 6, 4]]
    sections = [make_section(; v0 = d) for d in shuffled]
    sort!(sections; lt = (a, b) -> GO.edge_angle_compare(m, a, b; exact = True()) < 0)
    @test [GO.get_vertex(ns, 0) for ns in sections] == dirs_ccw

    # Equal-angle (collinear, same quadrant) sections compare equal.
    a = make_section(; v0 = (2.0, 2.0))
    b = make_section(; v0 = (1.0, 1.0))
    @test GO.edge_angle_compare(m, a, b; exact = True()) == 0

    # Symbolic crossing-node apex: comparator works without a constructed
    # intersection coordinate (design D2).
    xnode = GO.crossing_node((-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0))
    xdirs_ccw = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    xsections = [make_section(; node = xnode, v0 = d) for d in xdirs_ccw[[3, 1, 4, 2]]]
    sort!(xsections; lt = (a, b) -> GO.edge_angle_compare(m, a, b; exact = True()) < 0)
    @test [GO.get_vertex(ns, 0) for ns in xsections] == xdirs_ccw
end

@testset "compare_to" begin
    # A sorts before B
    @test GO.compare_to(make_section(; is_a = true), make_section(; is_a = false)) < 0
    @test GO.compare_to(make_section(; is_a = false), make_section(; is_a = true)) > 0
    # lines sort before areas
    @test GO.compare_to(make_section(; dim = GO.DIM_L), make_section(; dim = GO.DIM_A)) < 0
    # then id, then ring id
    @test GO.compare_to(make_section(; id = 1), make_section(; id = 2)) < 0
    @test GO.compare_to(make_section(; ring_id = 0), make_section(; ring_id = 1)) < 0
    # then edge vertices, `nothing` (Java null) sorting below non-null,
    # coordinates lexicographic on (x, y)
    @test GO.compare_to(make_section(; v0 = nothing), make_section(; v0 = (0.0, 0.0))) < 0
    @test GO.compare_to(make_section(; v0 = (0.0, 0.0)), make_section(; v0 = nothing)) > 0
    @test GO.compare_to(make_section(; v0 = (1.0, 0.0)), make_section(; v0 = (1.0, 2.0))) < 0
    @test GO.compare_to(make_section(; v0 = (1.0, 0.0)), make_section(; v0 = (2.0, 0.0))) < 0
    @test GO.compare_to(make_section(; v1 = nothing), make_section(; v1 = (0.0, 1.0))) < 0
    @test GO.compare_to(make_section(; v1 = (0.0, 1.0)), make_section(; v1 = (0.0, 2.0))) < 0
    @test GO.compare_to(make_section(; v1 = nothing), make_section(; v1 = nothing)) == 0
    @test GO.compare_to(make_section(), make_section()) == 0
    # geometry comparison dominates dimension, dimension dominates id, ...
    @test GO.compare_to(make_section(; is_a = true, dim = GO.DIM_A),
        make_section(; is_a = false, dim = GO.DIM_L)) < 0
    @test GO.compare_to(make_section(; dim = GO.DIM_L, id = 9),
        make_section(; dim = GO.DIM_A, id = 1)) < 0
    @test GO.compare_to(make_section(; id = 1, ring_id = 9),
        make_section(; id = 2, ring_id = 0)) < 0
end

@testset "NodeSections prepare_sections! ordering invariant" begin
    nss = GO.NodeSections(NODE)
    @test GO.get_coordinate(nss) === NODE
    s_line_a = make_section(; is_a = true, dim = GO.DIM_L, id = 5, ring_id = -1, v1 = nothing)
    s_shell_a1 = make_section(; is_a = true, dim = GO.DIM_A, id = 1, ring_id = 0)
    s_hole_a1 = make_section(; is_a = true, dim = GO.DIM_A, id = 1, ring_id = 1)
    s_shell_a2 = make_section(; is_a = true, dim = GO.DIM_A, id = 2, ring_id = 0)
    s_line_b = make_section(; is_a = false, dim = GO.DIM_L, id = 3, ring_id = -1, v1 = nothing)
    s_shell_b = make_section(; is_a = false, dim = GO.DIM_A, id = 1, ring_id = 0)
    for s in (s_shell_b, s_shell_a2, s_hole_a1, s_line_b, s_line_a, s_shell_a1)
        GO.add_node_section!(nss, s)
    end
    GO.prepare_sections!(nss)
    # A before B (dominating dimension: the B line sorts after every A area);
    # within a geometry lines before areas; same-polygon sections contiguous.
    @test nss.sections == [s_line_a, s_shell_a1, s_hole_a1, s_shell_a2, s_line_b, s_shell_b]
end

@testset "NodeSections" begin
    m = Planar()
    poly_a = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)]])
    poly_b = GI.Polygon([[(0.0, 0.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 0.0)]])

    nss = GO.NodeSections(NODE)
    @test !GO.has_interaction_ab(nss)
    # A-shell corner at the node (CW orientation: interior to the right)
    sa = make_section(; is_a = true, id = 1, poly = poly_a, v0 = (0.0, 1.0), v1 = (1.0, 0.0))
    GO.add_node_section!(nss, sa)
    @test !GO.has_interaction_ab(nss)
    # B-shell corner touching at the node from the opposite quadrant
    sb = make_section(; is_a = false, id = 1, poly = poly_b, v0 = (-1.0, 0.0), v1 = (0.0, -1.0))
    GO.add_node_section!(nss, sb)
    @test GO.has_interaction_ab(nss)

    @test GO.get_polygonal(nss, true) === poly_a
    @test GO.get_polygonal(nss, false) === poly_b

    # getPolygonal skips sections without a parent polygonal
    nss_line = GO.NodeSections(NODE)
    GO.add_node_section!(nss_line, make_section(; dim = GO.DIM_L, ring_id = -1, poly = nothing))
    @test GO.get_polygonal(nss_line, true) === nothing
    @test GO.get_polygonal(nss_line, false) === nothing
    @test !GO.has_interaction_ab(nss_line)

    # create_node on a simple two-area touch (full port as of Task 17:
    # returns the assembled RelateNode).
    #
    # Hand-trace: prepared order [sa, sb] (A before B). sa: enter (0,1)@90°
    # reverse {A: L=I,On=B,R=E}, exit (1,0)@0° forward {A: L=E,On=B,R=I};
    # the exit inserts BEFORE (0,1) in the CCW wheel (0° < 90°). No edges lie
    # strictly between, and the prev/next interior checks see EXTERIOR sides,
    # so no propagation. sb appends (-1,0)@180° {B: I,B,E} and (0,-1)@270°
    # {B: E,B,I}; the A-edge B-labels stay unknown until finish! (addEdges
    # never retro-labels other-geometry edges).
    node = GO.create_node(m, nss; exact = True())
    @test node isa GO.RelateNode
    edges = GO.get_edges(node)
    @test [e.dir_pt for e in edges] == [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    locs(e, isa_g) = (GO.location(e, isa_g, GO.POS_LEFT),
        GO.location(e, isa_g, GO.POS_ON), GO.location(e, isa_g, GO.POS_RIGHT))
    LI, LB, LE, LN = GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.LOC_NONE
    @test locs(edges[1], true) == (LE, LB, LI) && locs(edges[2], true) == (LI, LB, LE)
    @test locs(edges[3], false) == (LI, LB, LE) && locs(edges[4], false) == (LE, LB, LI)
    @test !GO.is_known(edges[1], false) && !GO.is_known(edges[3], true)
    @test GO.has_exterior_edge(node, true) && GO.has_exterior_edge(node, false)

    # Multiple sections of the same polygon route through PolygonNodeConverter
    # (Task 16): two shell corners of one polygon are rewritten to themselves,
    # in edge-angle order ((0,1) at 90° before (-1,0) at 180°).
    #
    # Hand-trace (Task 17): multi_1 (corner (0,1)→(1,0), interior sector CCW
    # 90°→0°, i.e. through 180° and 270°) builds wheel [(1,0){A:E,B,I},
    # (0,1){A:I,B,E}]. multi_2 (corner (-1,0)→(0,-1), sector 180°→270°)
    # appends (-1,0)@180° {I,B,E} and (0,-1)@270° {E,B,I}; then
    # updateIfAreaPrev(enter=3): prev edge (0,1) has LEFT=I → edge (-1,0)
    # becomes all-INTERIOR; updateIfAreaNext(exit=4): next edge (wraps to
    # (1,0)) has RIGHT=I → edge (0,-1) becomes all-INTERIOR. That matches the
    # geometry: multi_2's corner lies inside multi_1's interior sector.
    nss_multi = GO.NodeSections(NODE)
    multi_1 = make_section(; id = 1, v0 = (0.0, 1.0), v1 = (1.0, 0.0))
    multi_2 = make_section(; id = 1, v0 = (-1.0, 0.0), v1 = (0.0, -1.0))
    GO.add_node_section!(nss_multi, multi_1)
    GO.add_node_section!(nss_multi, multi_2)
    node_multi = GO.create_node(m, nss_multi; exact = True())
    edges_multi = GO.get_edges(node_multi)
    @test [e.dir_pt for e in edges_multi] ==
        [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
    @test locs(edges_multi[1], true) == (LE, LB, LI)
    @test locs(edges_multi[2], true) == (LI, LB, LE)
    @test locs(edges_multi[3], true) == (LI, LI, LI)
    @test locs(edges_multi[4], true) == (LI, LI, LI)
end

# Port of JTS PolygonNodeConverterTest.java — every test method, plus the
# checkConversion / checkSectionsEqual / sort / section helpers. The Java
# section helpers build A-geometry area sections of polygon 1 at node (5,5)
# (no parent polygonal, node not at a vertex); equality of section lists is
# up to edge-angle order, compared with NodeSection.compareTo.
@testset "PolygonNodeConverter" begin
    m = Planar()

    # Port of section(ringId, v0x, v0y, nx, ny, v1x, v1y).
    section(ring_id, v0x, v0y, nx, ny, v1x, v1y) = GO.NodeSection(
        true, GO.DIM_A, Int32(1), Int32(ring_id), nothing, false,
        (Float64(v0x), Float64(v0y)),
        GO.vertex_node((Float64(nx), Float64(ny))),
        (Float64(v1x), Float64(v1y)))
    # Ports of sectionShell / sectionHole.
    section_shell(coords...) = section(0, coords...)
    section_hole(coords...) = section(1, coords...)

    # Port of sort(List<NodeSection>): EdgeAngleComparator ordering.
    sort_sections!(ns) =
        sort!(ns; lt = (a, b) -> GO.edge_angle_compare(m, a, b; exact = True()) < 0)

    # Port of checkSectionsEqual.
    function is_sections_equal(ns1, ns2)
        length(ns1) == length(ns2) || return false
        sort_sections!(ns1)
        sort_sections!(ns2)
        for i in eachindex(ns1)
            GO.compare_to(ns1[i], ns2[i]) == 0 || return false
        end
        return true
    end

    # Port of checkConversion.
    function check_conversion(input, expected)
        actual = GO.polygon_node_convert(m, input; exact = True())
        @test is_sections_equal(actual, expected)
    end

    @testset "testShells" begin
        check_conversion(
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 9, 9),
                section_shell(8, 9, 5, 5, 6, 9),
                section_shell(4, 9, 5, 5, 2, 9)],
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 9, 9),
                section_shell(8, 9, 5, 5, 6, 9),
                section_shell(4, 9, 5, 5, 2, 9)])
    end

    @testset "testShellAndHole" begin
        check_conversion(
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 9, 9),
                section_hole(6, 0, 5, 5, 4, 0)],
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 4, 0),
                section_shell(6, 0, 5, 5, 9, 9)])
    end

    @testset "testShellsAndHoles" begin
        check_conversion(
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 9, 9),
                section_hole(6, 0, 5, 5, 4, 0),
                section_shell(8, 8, 5, 5, 1, 8),
                section_hole(4, 8, 5, 5, 6, 8)],
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 4, 0),
                section_shell(6, 0, 5, 5, 9, 9),
                section_shell(4, 8, 5, 5, 1, 8),
                section_shell(8, 8, 5, 5, 6, 8)])
    end

    @testset "testShellAnd2Holes" begin
        check_conversion(
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 9, 9),
                section_hole(7, 0, 5, 5, 6, 0),
                section_hole(4, 0, 5, 5, 3, 0)],
            GO.NodeSection[
                section_shell(1, 1, 5, 5, 3, 0),
                section_shell(4, 0, 5, 5, 6, 0),
                section_shell(7, 0, 5, 5, 9, 9)])
    end

    @testset "testHoles" begin
        check_conversion(
            GO.NodeSection[
                section_hole(7, 0, 5, 5, 6, 0),
                section_hole(4, 0, 5, 5, 3, 0)],
            GO.NodeSection[
                section_shell(4, 0, 5, 5, 6, 0),
                section_shell(7, 0, 5, 5, 3, 0)])
    end
end

# Tests for relate_node.jl (Task 17): ports of JTS RelateEdge.java and
# RelateNode.java. No dedicated JUnit file exists; each configuration below
# was hand-traced against the Java semantics (trace reasoning in comments).
#
# Conventions used in the traces (derived from the Java sources):
# - Area sections arrive in canonical orientation: CW shells / CCW holes,
#   i.e. polygon interior on the RIGHT of travel v0 → node → v1. The interior
#   sector at the node therefore spans CCW from ray node→v0 to ray node→v1.
# - RelateNode.addEdges(A-section): the entering edge (dirPt v0) is added
#   with isForward=false → LEFT=INTERIOR, RIGHT=EXTERIOR; the exiting edge
#   (dirPt v1) with isForward=true → LEFT=EXTERIOR, RIGHT=INTERIOR; both get
#   ON=BOUNDARY. (RelateEdge.setLocationsArea.)
# - Line edges get LEFT=RIGHT=EXTERIOR, ON=INTERIOR, dim L.
#   (RelateEdge.setLocationsLine.)
# - The wheel is kept sorted CCW by angle from the positive X-axis
#   (RelateNode.addEdge insertion via compareToEdge).
@testset "RelateEdge + RelateNode" begin
    m = Planar()
    LI, LB, LE, LN = GO.LOC_INTERIOR, GO.LOC_BOUNDARY, GO.LOC_EXTERIOR, GO.LOC_NONE
    locs(e, isa_g) = (GO.location(e, isa_g, GO.POS_LEFT),
        GO.location(e, isa_g, GO.POS_ON), GO.location(e, isa_g, GO.POS_RIGHT))
    dirs(node) = [e.dir_pt for e in GO.get_edges(node)]

    # Build a node through the full NodeSections.createNode pipeline.
    function build_node(sections...; node_key = NODE)
        nss = GO.NodeSections(node_key)
        for s in sections
            GO.add_node_section!(nss, s)
        end
        return GO.create_node(m, nss; exact = True())
    end

    line_section(; kw...) = make_section(; dim = GO.DIM_L, ring_id = -1, kw...)

    @testset "RelateEdge basics" begin
        # Factory: area dim → area edge (sides labeled by direction);
        # any other dim → line edge.
        ea = GO.relate_edge(NODE, (1.0, 0.0), true, GO.DIM_A, true)
        @test ea.a_dim == GO.DIM_A
        @test locs(ea, true) == (LE, LB, LI)        # forward: interior on R
        @test locs(ea, false) == (LN, LN, LN)       # B untouched
        @test GO.is_known(ea, true) && !GO.is_known(ea, false)

        er = GO.relate_edge(NODE, (1.0, 0.0), true, GO.DIM_A, false)
        @test locs(er, true) == (LI, LB, LE)        # reverse: interior on L

        el = GO.relate_edge(NODE, (1.0, 0.0), false, GO.DIM_L, true)
        @test el.b_dim == GO.DIM_L && el.a_dim == GO.DIM_UNKNOWN_EDGE
        @test locs(el, false) == (LE, LI, LE)

        # Explicit-locations constructor (Java RelateEdge(node, pt, isA,
        # locLeft, locRight, locLine)) forces dim 2.
        ex = GO.RelateEdge(NODE, (1.0, 0.0), true, GO.LOC_INTERIOR,
            GO.LOC_INTERIOR, GO.LOC_INTERIOR)
        @test ex.a_dim == GO.DIM_A && locs(ex, true) == (LI, LI, LI)

        # location / is_interior / set_location! / set_all_locations! /
        # set_unknown_locations!
        @test GO.is_interior(ea, true, GO.POS_RIGHT)
        @test !GO.is_interior(ea, true, GO.POS_LEFT)
        @test_throws ArgumentError GO.location(ea, true, Int8(99))
        GO.set_location!(ea, false, GO.POS_ON, GO.LOC_EXTERIOR)
        @test GO.location(ea, false, GO.POS_ON) == LE
        GO.set_unknown_locations!(ea, false, GO.LOC_INTERIOR)
        @test locs(ea, false) == (LI, LE, LI)       # ON was known, kept
        GO.set_all_locations!(ea, false, GO.LOC_EXTERIOR)
        @test locs(ea, false) == (LE, LE, LE)
        GO.set_dim_locations!(ea, false, GO.DIM_L, GO.LOC_INTERIOR)
        @test ea.b_dim == GO.DIM_L && locs(ea, false) == (LI, LI, LI)
        GO.set_area_interior!(er, true)
        @test locs(er, true) == (LI, LI, LI)
        @test er.a_dim == GO.DIM_A                  # dim untouched

        # statics: find_known_edge_index (1-based, 0 if none) and
        # set_all_area_interior!
        edges = [GO.relate_edge(NODE, (1.0, 0.0), false, GO.DIM_L, true),
            GO.relate_edge(NODE, (0.0, 1.0), true, GO.DIM_A, true)]
        @test GO.find_known_edge_index(edges, true) == 2
        @test GO.find_known_edge_index(edges, false) == 1
        @test GO.find_known_edge_index([edges[2]], false) == 0
        GO.set_all_area_interior!(edges, true)
        @test locs(edges[1], true) == (LI, LI, LI)
        @test locs(edges[2], true) == (LI, LI, LI)

        # compare_to_edge: CCW angle comparison around the node
        # (negative: this edge below the query direction; positive: above).
        e0 = GO.relate_edge(NODE, (1.0, 1.0), true, GO.DIM_A, true)  # 45°
        @test GO.compare_to_edge(m, e0, (0.0, 1.0); exact = True()) < 0   # 45 < 90
        @test GO.compare_to_edge(m, e0, (1.0, 0.0); exact = True()) > 0   # 45 > 0
        @test GO.compare_to_edge(m, e0, (2.0, 2.0); exact = True()) == 0  # collinear
        # toString port (smoke)
        @test sprint(show, e0) isa String
    end

    # Configuration 1: two lines crossing at a shared vertex (0,0).
    #
    # Hand-trace: sections (sorted A first) are A-line v0=(-1,0), v1=(1,0)
    # and B-line v0=(0,-1), v1=(0,1); dim L → addLineEdge per vertex.
    # A: (-1,0)@180° starts the wheel; (1,0)@0° compares "further" (+1
    # against 180°) and inserts before it → [(1,0), (-1,0)].
    # B: (0,-1)@270° appends; (0,1)@90° inserts before (-1,0) →
    # [(1,0)@0°, (0,1)@90°, (-1,0)@180°, (0,-1)@270°].
    # Every edge: own geometry dim L with (E, I, E), other geometry fully
    # unknown — no area labels anywhere. finish!(node, false, false) then
    # propagates: each edge's unknown side/on locations become the CCW-prior
    # known LEFT location, which is EXTERIOR everywhere here.
    @testset "config 1: crossing lines at a vertex node" begin
        sA = line_section(; v0 = (-1.0, 0.0), v1 = (1.0, 0.0))
        sB = line_section(; is_a = false, v0 = (0.0, -1.0), v1 = (0.0, 1.0))
        node = build_node(sA, sB)
        edges = GO.get_edges(node)
        @test length(edges) == 4
        @test dirs(node) == [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
        # all dims L for the owning geometry, unknown for the other
        @test edges[1].a_dim == GO.DIM_L && edges[3].a_dim == GO.DIM_L
        @test edges[2].b_dim == GO.DIM_L && edges[4].b_dim == GO.DIM_L
        @test edges[1].b_dim == GO.DIM_UNKNOWN_EDGE
        @test edges[2].a_dim == GO.DIM_UNKNOWN_EDGE
        # line labels: (E, I, E) for own geometry; nothing known for other
        for (i, e) in pairs(edges)
            own = isodd(i)   # edges 1,3 are A's; 2,4 are B's
            @test locs(e, own) == (LE, LI, LE)
            @test locs(e, !own) == (LN, LN, LN)
        end
        # no area labels: no BOUNDARY/INTERIOR side locations at all
        @test all(e -> GO.location(e, true, GO.POS_LEFT) != LB &&
            GO.location(e, false, GO.POS_LEFT) != LB, edges)

        GO.finish!(node, false, false)
        for (i, e) in pairs(edges)
            own = isodd(i)
            @test locs(e, own) == (LE, LI, LE)      # known labels survive
            @test locs(e, !own) == (LE, LE, LE)     # unknowns → EXTERIOR
        end
        @test GO.has_exterior_edge(node, true) && GO.has_exterior_edge(node, false)

        # Same configuration at a *symbolic* crossing node (design D2): the
        # wheel is ordered around the crossing of (-1,0)-(1,0) × (0,-1)-(0,1)
        # without any constructed apex coordinate. (The zero-length-edge
        # guard must not engage for crossing keys.)
        xnode = GO.crossing_node((-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0))
        sAx = line_section(; at_vertex = false, v0 = (-1.0, 0.0), v1 = (1.0, 0.0), node = xnode)
        sBx = line_section(; is_a = false, at_vertex = false,
            v0 = (0.0, -1.0), v1 = (0.0, 1.0), node = xnode)
        xn = build_node(sAx, sBx; node_key = xnode)
        @test dirs(xn) == [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)]
        @test locs(GO.get_edges(xn)[1], true) == (LE, LI, LE)
        @test locs(GO.get_edges(xn)[2], false) == (LE, LI, LE)
    end

    # Configuration 2: an area corner — two edges of one CW-shell polygon.
    #
    # Hand-trace: corner v0=(1,0), v1=(0,1) (CW shell through the node, e.g.
    # ring … → (1,0) → (0,0) → (0,1) → …; interior is the first quadrant,
    # the sector CCW from ray node→v0 @0° to ray node→v1 @90°).
    # addEdges: entering edge (1,0) reverse → {A: L=I, On=B, R=E}; exiting
    # edge (0,1) forward appends after 0° → {A: L=E, On=B, R=I}.
    # updateEdgesInArea(1→2) walks no edges (adjacent); prev/next interior
    # checks see EXTERIOR → no propagation. Interior ends up on the correct
    # side: LEFT of the entering ray (CCW side toward the sector) and RIGHT
    # of the exiting ray.
    @testset "config 2: area corner of one polygon" begin
        corner = make_section(; v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        node = build_node(corner)
        edges = GO.get_edges(node)
        @test length(edges) == 2
        @test dirs(node) == [(1.0, 0.0), (0.0, 1.0)]
        @test edges[1].a_dim == GO.DIM_A && edges[2].a_dim == GO.DIM_A
        @test locs(edges[1], true) == (LI, LB, LE)
        @test locs(edges[2], true) == (LE, LB, LI)
        @test locs(edges[1], false) == (LN, LN, LN)
        @test GO.has_exterior_edge(node, true)
        @test !GO.has_exterior_edge(node, false)    # nothing known for B
    end

    # Configuration 3: area corner of A + line end of B at the same node.
    #
    # Hand-trace: A corner as config 2 (interior sector 0°→90°). B line ends
    # at the node arriving from (1,1): section v0=(1,1)@45°, v1=nothing.
    # Sorted A before B, so the wheel is [(1,0), (0,1)] before the B line
    # edge inserts between them (0° < 45° < 90°) → [(1,0), (1,1), (0,1)];
    # the nothing vertex is skipped (Java addEdge null guard). Pre-finish
    # the line edge knows nothing about A. finish!(node, false, false):
    # propagateSideLocations(A) starts at edge 1 (first A-known), carries
    # currLoc = LEFT(A) of (1,0) = INTERIOR onto the line edge → its A
    # locations all become INTERIOR (the line end lies inside A); then
    # currLoc resets to EXTERIOR past (0,1)'s LEFT... (0,1) is fully known,
    # so nothing else changes. propagateSideLocations(B) starts at the line
    # edge, carries its LEFT(B) = EXTERIOR onto both area edges.
    @testset "config 3: area corner of A + line end of B" begin
        corner = make_section(; v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        line_end = line_section(; is_a = false, v0 = (1.0, 1.0), v1 = nothing)
        node = build_node(corner, line_end)
        edges = GO.get_edges(node)
        @test length(edges) == 3
        @test dirs(node) == [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        line_edge = edges[2]
        @test line_edge.b_dim == GO.DIM_L
        @test locs(line_edge, false) == (LE, LI, LE)
        @test locs(line_edge, true) == (LN, LN, LN)     # A unknown pre-finish

        GO.finish!(node, false, false)
        # the line edge is interior to A on every position (line end is
        # inside A's interior sector); its dimension for A stays unknown
        @test locs(line_edge, true) == (LI, LI, LI)
        @test line_edge.a_dim == GO.DIM_UNKNOWN_EDGE
        # the area edges are exterior to B everywhere
        @test locs(edges[1], false) == (LE, LE, LE)
        @test locs(edges[3], false) == (LE, LE, LE)
        # A labels on its own edges survive finish!
        @test locs(edges[1], true) == (LI, LB, LE)
        @test locs(edges[3], true) == (LE, LB, LI)

        # Same-geometry variant: when the line and the area belong to ONE
        # geometry, the labeling happens during addEdges (updateEdgesInArea)
        # rather than at finish!. prepareSections puts the line first, so
        # the wheel is [(1,1) line] when the corner is added; the corner's
        # entering edge (1,0) inserts at 1, exiting edge (0,1) appends →
        # [(1,0), (1,1), (0,1)] with index0=1, index1=3, and
        # updateEdgesInArea marks edge 2 (the line edge, strictly between
        # the corner edges in CCW order) all-INTERIOR for A. Its dim stays L
        # (setAreaInterior only touches locations).
        line_a = line_section(; v0 = (1.0, 1.0), v1 = nothing, id = 2)
        node2 = build_node(corner, line_a)
        edges2 = GO.get_edges(node2)
        @test dirs(node2) == [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        @test locs(edges2[2], true) == (LI, LI, LI)
        @test edges2[2].a_dim == GO.DIM_L
    end

    # Configuration 4: two area corners (A and B) with coincident edges —
    # the collinear-edge merge case, plus the area-over-line dim override.
    #
    # Hand-trace (coincident A/B corners, both v0=(1,0), v1=(0,1)):
    # A corner builds [(1,0){A: I,B,E}, (0,1){A: E,B,I}]. B corner's
    # entering edge compares equal (comp == 0) to (1,0) → RelateEdge.merge
    # with isA=false: B is unknown on that edge, so the !isKnown branch sets
    # B dim/locations directly: {B: I,B,E}. Likewise the exiting edge merges
    # into (0,1) → {B: E,B,I}. The wheel stays at 2 edges, both fully
    # labeled for A and B identically (the boundaries coincide).
    @testset "config 4: coincident area corners of A and B" begin
        corner_a = make_section(; v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        corner_b = make_section(; is_a = false, v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        node = build_node(corner_a, corner_b)
        edges = GO.get_edges(node)
        @test length(edges) == 2
        @test dirs(node) == [(1.0, 0.0), (0.0, 1.0)]
        @test edges[1].a_dim == GO.DIM_A && edges[1].b_dim == GO.DIM_A
        @test locs(edges[1], true) == (LI, LB, LE)
        @test locs(edges[1], false) == (LI, LB, LE)
        @test locs(edges[2], true) == (LE, LB, LI)
        @test locs(edges[2], false) == (LE, LB, LI)

        # Area-over-line dim override (RelateEdge.mergeDimEdgeLoc): geometry
        # B contributes a line edge and a coincident area edge. The line
        # sorts first (prepareSections: lines before areas), so the wheel
        # holds (1,0){B: dim L, E,I,E} when the corner's entering edge
        # merges into it: isKnown(B) → mergeDimEdgeLoc upgrades dim L → A
        # and ON → BOUNDARY; mergeSideLocation sets LEFT to INTERIOR (curr
        # EXTERIOR yields) and leaves RIGHT EXTERIOR.
        line_b = line_section(; is_a = false, id = 2, v0 = (1.0, 0.0), v1 = nothing)
        corner_b2 = make_section(; is_a = false, v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        node2 = build_node(line_b, corner_b2)
        edges2 = GO.get_edges(node2)
        @test length(edges2) == 2
        merged = edges2[1]
        @test merged.dir_pt == (1.0, 0.0)
        @test merged.b_dim == GO.DIM_A              # dim override L → A
        @test locs(merged, false) == (LI, LB, LE)   # ON overridden to BOUNDARY
        @test locs(edges2[2], false) == (LE, LB, LI)

        # Merge side-location precedence: INTERIOR wins over EXTERIOR. Two
        # self-touch corners of one polygon sharing the ray (0,1): corner 1
        # (1,0)→(0,1) (sector 0°→90°) and corner 2 (0,1)→(-1,0) (sector
        # 90°→180°). After the converter (pass-through for two shells) the
        # corners are added in v0-angle order: corner 1 builds
        # [(1,0){I,B,E}, (0,1){E,B,I}]; corner 2's entering edge merges into
        # (0,1) reverse → LEFT: EXTERIOR→INTERIOR, RIGHT: INTERIOR kept →
        # {I,B,I}; the exiting edge (-1,0) appends {E,B,I}. Then
        # updateIfAreaPrev(enter=2): prev edge (1,0) has LEFT=I → the shared
        # edge (0,1) becomes all-INTERIOR — the coincident boundary ray is
        # interior to the union of the two corners (0°→180°), exactly the
        # self-touch (maximal ring) semantics.
        corner_1 = make_section(; v0 = (1.0, 0.0), v1 = (0.0, 1.0))
        corner_2 = make_section(; v0 = (0.0, 1.0), v1 = (-1.0, 0.0))
        node3 = build_node(corner_1, corner_2)
        edges3 = GO.get_edges(node3)
        @test length(edges3) == 3
        @test dirs(node3) == [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)]
        @test locs(edges3[1], true) == (LI, LB, LE)
        @test locs(edges3[2], true) == (LI, LI, LI)     # shared ray: interior
        @test locs(edges3[3], true) == (LE, LB, LI)
    end

    # Carried-over from the Task 16 review: a node mixing a converted
    # multi-section polygon group (shell + hole of A's polygon 1) with
    # another geometry's singleton corner (B).
    #
    # Hand-trace: A shell runs straight through the node, v0=(0,-1)@270°,
    # v1=(0,1)@90° (interior = east half-plane). A hole touches the node
    # occupying the thin east wedge 315°→45°: CCW hole corner v0=(1,1)@45°,
    # v1=(1,-1)@315° (polygon interior CCW from 45° to 315° — the rule is
    # the same right-of-travel rule as for shells). PolygonNodeConverter
    # (sorted by v0 angle: hole@45°, shell@270°; findShell → shell) rewrites
    # them into the self-touch corners (shell.v0 → hole.v1) = (0,-1)→(1,-1)
    # [sector 270°→315°] and (hole.v0 → shell.v1) = (1,1)→(0,1) [sector
    # 45°→90°] — partitioning the polygon interior near the node. B's
    # CW corner v0=(-1,1)@135°, v1=(-1,-1)@225° (interior = the west wedge,
    # inside A's exterior) is added as a singleton, untouched by the
    # converter. Wheel insertion order gives CCW dirs
    # [(1,1)@45°, (0,1)@90°, (-1,1)@135°, (-1,-1)@225°, (0,-1)@270°,
    # (1,-1)@315°]; no updateEdgesInArea/prev/next propagation fires (each
    # corner's edges are adjacent in the wheel, and neighboring edges of the
    # other corners show EXTERIOR or unknown sides at the decisive checks).
    # finish!(node, false, false) then labels A's edges exterior-for-B and
    # B's edges exterior-for-A (each lies in the other's exterior).
    @testset "converted group + other polygon's singleton" begin
        shell_a = make_section(; id = 1, ring_id = 0, v0 = (0.0, -1.0), v1 = (0.0, 1.0))
        hole_a = make_section(; id = 1, ring_id = 1, v0 = (1.0, 1.0), v1 = (1.0, -1.0))
        corner_b = make_section(; is_a = false, id = 1, v0 = (-1.0, 1.0), v1 = (-1.0, -1.0))
        node = build_node(shell_a, hole_a, corner_b)
        edges = GO.get_edges(node)
        @test length(edges) == 6
        @test dirs(node) == [(1.0, 1.0), (0.0, 1.0), (-1.0, 1.0),
            (-1.0, -1.0), (0.0, -1.0), (1.0, -1.0)]
        # A's converted self-touch corners: (0,-1)→(1,-1) and (1,1)→(0,1)
        @test locs(edges[1], true) == (LI, LB, LE)      # (1,1): enter of corner 2
        @test locs(edges[2], true) == (LE, LB, LI)      # (0,1): exit of corner 2
        @test locs(edges[5], true) == (LI, LB, LE)      # (0,-1): enter of corner 1
        @test locs(edges[6], true) == (LE, LB, LI)      # (1,-1): exit of corner 1
        # B's singleton corner
        @test locs(edges[3], false) == (LI, LB, LE)     # (-1,1): enter
        @test locs(edges[4], false) == (LE, LB, LI)     # (-1,-1): exit
        # cross-geometry labels unknown before finish!
        @test locs(edges[1], false) == (LN, LN, LN)
        @test locs(edges[3], true) == (LN, LN, LN)

        GO.finish!(node, false, false)
        # B's corner lies in A's exterior and vice versa
        @test locs(edges[3], true) == (LE, LE, LE)
        @test locs(edges[4], true) == (LE, LE, LE)
        @test locs(edges[1], false) == (LE, LE, LE)
        @test locs(edges[5], false) == (LE, LE, LE)
        @test GO.has_exterior_edge(node, true) && GO.has_exterior_edge(node, false)
    end

    # finish! with isAreaInterior: in a mixed GC the node may lie in the
    # interior of an area of the other geometry — every edge then becomes
    # all-INTERIOR for that geometry (RelateNode.finishNode true branch).
    @testset "finish! isAreaInterior" begin
        sA = line_section(; v0 = (-1.0, 0.0), v1 = (1.0, 0.0))
        sB = line_section(; is_a = false, v0 = (0.0, -1.0), v1 = (0.0, 1.0))
        node = build_node(sA, sB)
        GO.finish!(node, true, false)
        for e in GO.get_edges(node)
            @test locs(e, true) == (LI, LI, LI)
        end
        @test !GO.has_exterior_edge(node, true)
        @test GO.has_exterior_edge(node, false)
    end
end

@testset "no GO methods on Base.merge!" begin
    # RelateEdge label merging must not extend Base.merge! (collection
    # semantics). It lives on the internal function `merge_edge!` instead.
    @test !any(m -> m.module == GO, methods(Base.merge!))
end
