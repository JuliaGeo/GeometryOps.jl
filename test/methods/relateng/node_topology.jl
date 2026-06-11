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

    # create_node on a simple two-area touch. Partial port (see
    # node_sections.jl): RelateNode lands in Task 17, so create_node returns
    # the prepared, converted section list the node's addEdges will consume —
    # one section per polygon here, A first. (Edge-count and label assertions
    # land with RelateNode in Task 17.)
    out = GO.create_node(m, nss; exact = True())
    @test length(out) == 2
    @test out[1] === sa
    @test out[2] === sb

    # Multiple sections of the same polygon route through PolygonNodeConverter
    # (Task 16): two shell corners of one polygon are rewritten to themselves,
    # in edge-angle order ((0,1) at 90° before (-1,0) at 180°).
    nss_multi = GO.NodeSections(NODE)
    multi_1 = make_section(; id = 1, v0 = (0.0, 1.0), v1 = (1.0, 0.0))
    multi_2 = make_section(; id = 1, v0 = (-1.0, 0.0), v1 = (0.0, -1.0))
    GO.add_node_section!(nss_multi, multi_1)
    GO.add_node_section!(nss_multi, multi_2)
    out_multi = GO.create_node(m, nss_multi; exact = True())
    @test length(out_multi) == 2
    @test GO.compare_to(out_multi[1], multi_1) == 0
    @test GO.compare_to(out_multi[2], multi_2) == 0
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
