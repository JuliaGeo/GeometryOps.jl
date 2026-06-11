# Tests for the RelateNG input facade (relate_geometry.jl). The
# "RelateGeometry" testset ports JTS RelateGeometryTest.java in full. JTS has
# no RelateSegmentStringTest, so the "RelateSegmentString" testset is
# hand-written against the Java RelateSegmentString.java semantics.

using Test
import GeometryOps as GO
import GeometryOps: Planar, True
import GeoInterface as GI
import LibGEOS as LG  # only for POLYGON EMPTY — GI wrappers cannot be empty
import Extents

relate_geom(geom) = GO.RelateGeometry(Planar(), geom; exact = True())

# Port of RelateGeometryTest.checkDimension.
function check_dimension(geom, expected_dim, expected_dim_real)
    rgeom = relate_geom(geom)
    @test GO.get_dimension(rgeom) == expected_dim
    @test GO.get_dimension_real(rgeom) == expected_dim_real
end

@testset "RelateGeometry" begin
    # testUniquePoints: MULTIPOINT ((0 0), (5 5), (5 0), (0 0))
    @testset "unique points" begin
        geom = GI.MultiPoint([(0.0, 0.0), (5.0, 5.0), (5.0, 0.0), (0.0, 0.0)])
        rgeom = relate_geom(geom)
        pts = GO.get_unique_points(rgeom)
        @test length(pts) == 3  # "Unique pts size"
    end

    # testBoundary: MULTILINESTRING ((0 0, 9 9), (9 9, 5 1))
    @testset "boundary" begin
        geom = GI.MultiLineString([
            GI.LineString([(0.0, 0.0), (9.0, 9.0)]),
            GI.LineString([(9.0, 9.0), (5.0, 1.0)]),
        ])
        rgeom = relate_geom(geom)
        @test GO.has_boundary(rgeom)  # "hasBoundary"
    end

    # testHasDimension:
    # GEOMETRYCOLLECTION (POLYGON ((1 9, 5 9, 5 5, 1 5, 1 9)),
    #                     LINESTRING (1 1, 5 4), POINT (6 5))
    @testset "has dimension" begin
        geom = GI.GeometryCollection([
            GI.Polygon([[(1.0, 9.0), (5.0, 9.0), (5.0, 5.0), (1.0, 5.0), (1.0, 9.0)]]),
            GI.LineString([(1.0, 1.0), (5.0, 4.0)]),
            GI.Point(6.0, 5.0),
        ])
        rgeom = relate_geom(geom)
        @test GO.has_dimension(rgeom, 0)  # "hasDimension 0"
        @test GO.has_dimension(rgeom, 1)  # "hasDimension 1"
        @test GO.has_dimension(rgeom, 2)  # "hasDimension 2"
    end

    # testDimension
    @testset "dimension" begin
        # POINT (0 0)
        check_dimension(GI.Point(0.0, 0.0), 0, 0)
        # LINESTRING (0 0, 0 0) — zero-length line is effectively a point
        check_dimension(GI.LineString([(0.0, 0.0), (0.0, 0.0)]), 1, 0)
        # LINESTRING (0 0, 9 9)
        check_dimension(GI.LineString([(0.0, 0.0), (9.0, 9.0)]), 1, 1)
        # LINESTRING (0 0, 0 0, 9 9)
        check_dimension(GI.LineString([(0.0, 0.0), (0.0, 0.0), (9.0, 9.0)]), 1, 1)
        # POLYGON ((1 9, 5 9, 5 5, 1 5, 1 9))
        check_dimension(
            GI.Polygon([[(1.0, 9.0), (5.0, 9.0), (5.0, 5.0), (1.0, 5.0), (1.0, 9.0)]]),
            2, 2)
        # GEOMETRYCOLLECTION (POLYGON ((1 9, 5 9, 5 5, 1 5, 1 9)),
        #                     LINESTRING (1 1, 5 4), POINT (6 5))
        check_dimension(GI.GeometryCollection([
            GI.Polygon([[(1.0, 9.0), (5.0, 9.0), (5.0, 5.0), (1.0, 5.0), (1.0, 9.0)]]),
            GI.LineString([(1.0, 1.0), (5.0, 4.0)]),
            GI.Point(6.0, 5.0),
        ]), 2, 2)
        # GEOMETRYCOLLECTION (POLYGON EMPTY, LINESTRING (1 1, 5 4), POINT (6 5))
        # — the empty polygon still counts for getDimension (Java semantics),
        # but not for the real (non-empty) dimension.
        check_dimension(GI.GeometryCollection([
            LG.readgeom("POLYGON EMPTY"),
            GI.LineString([(1.0, 1.0), (5.0, 4.0)]),
            GI.Point(6.0, 5.0),
        ]), 2, 1)
    end
end

# The remaining testsets have no JTS JUnit counterpart; each case is
# hand-verified against RelateGeometry.java / RelateSegmentString.java.

@testset "RelateGeometry predicates" begin
    poly = GI.Polygon([[(1.0, 9.0), (5.0, 9.0), (5.0, 5.0), (1.0, 5.0), (1.0, 9.0)]])
    line = GI.LineString([(1.0, 1.0), (5.0, 4.0)])
    pt = GI.Point(6.0, 5.0)

    @testset "is_polygonal" begin
        @test GO.is_polygonal(relate_geom(poly))
        @test GO.is_polygonal(relate_geom(GI.MultiPolygon([poly])))
        @test !GO.is_polygonal(relate_geom(GI.GeometryCollection([poly])))
        @test !GO.is_polygonal(relate_geom(line))
    end

    @testset "has_edges / has_area_and_line" begin
        @test GO.has_edges(relate_geom(poly))
        @test GO.has_edges(relate_geom(line))
        @test !GO.has_edges(relate_geom(pt))
        @test GO.has_area_and_line(relate_geom(GI.GeometryCollection([poly, line])))
        @test !GO.has_area_and_line(relate_geom(poly))
    end

    @testset "is_self_noding_required" begin
        # points and polygonal geometries never require self-noding
        @test !GO.is_self_noding_required(relate_geom(pt))
        @test !GO.is_self_noding_required(relate_geom(GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)])))
        @test !GO.is_self_noding_required(relate_geom(poly))
        @test !GO.is_self_noding_required(relate_geom(GI.MultiPolygon([poly])))
        # lines may self-cross
        @test GO.is_self_noding_required(relate_geom(line))
        # a GC with a single polygon does not need noding
        @test !GO.is_self_noding_required(relate_geom(GI.GeometryCollection([poly])))
        # GCs with only points do not need noding
        @test !GO.is_self_noding_required(relate_geom(
            GI.GeometryCollection([pt, GI.Point(7.0, 5.0)])))
        # mixed-dimension GCs may have overlapping elements
        @test GO.is_self_noding_required(relate_geom(GI.GeometryCollection([poly, line])))
    end

    @testset "get_effective_points" begin
        # Points covered by a higher-dimension element are not effective.
        gc = GI.GeometryCollection([
            GI.Point(5.0, 5.0),    # inside the polygon — covered
            GI.Point(20.0, 20.0),  # outside — effective
            GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]]),
        ])
        rgeom = relate_geom(gc)
        eff = GO.get_effective_points(rgeom)
        @test length(eff) == 1
        @test GO._tuple_point(only(eff)) == (20.0, 20.0)
        # For a P-dimension geometry all points are returned unfiltered.
        mp = GI.MultiPoint([(0.0, 0.0), (5.0, 5.0)])
        @test length(GO.get_effective_points(relate_geom(mp))) == 2
    end
end

@testset "RelateSegmentString" begin
    # POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (2 2, 2 4, 4 4, 4 2, 2 2))
    # shell is CCW (must be flipped to CW), hole is CW (must be flipped to CCW)
    poly = GI.Polygon([
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
        [(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0), (2.0, 2.0)],
    ])

    @testset "polygon ring extraction" begin
        rgeom = relate_geom(poly)
        sss = GO.extract_segment_strings(rgeom, GO.GEOM_A, nothing)
        @test length(sss) == 2
        shell, hole = sss
        @test shell.is_a && hole.is_a
        @test shell.dim == GO.DIM_A && hole.dim == GO.DIM_A
        @test shell.id == 1 && hole.id == 1
        @test shell.ring_id == 0 && hole.ring_id == 1
        @test shell.parent_polygonal === poly && hole.parent_polygonal === poly
        @test shell.input_geom === rgeom
        # shell reoriented CW, hole reoriented CCW
        @test shell.pts[1] == (0.0, 0.0) && shell.pts[2] == (0.0, 10.0)
        @test hole.pts[1] == (2.0, 2.0) && hole.pts[2] == (4.0, 2.0)
        @test GO.is_closed(shell) && GO.is_closed(hole)
    end

    @testset "multipolygon parent and element ids" begin
        polys = [
            GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]),
            GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 5.0)]]),
        ]
        mp = GI.MultiPolygon(polys)
        rgeom = relate_geom(mp)
        sss = GO.extract_segment_strings(rgeom, GO.GEOM_B, nothing)
        @test length(sss) == 2
        @test !sss[1].is_a && !sss[2].is_a
        @test sss[1].id == 1 && sss[2].id == 2
        @test sss[1].parent_polygonal === mp && sss[2].parent_polygonal === mp
    end

    @testset "extent filter" begin
        polys = [
            GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]),
            GI.Polygon([[(50.0, 50.0), (60.0, 50.0), (60.0, 60.0), (50.0, 50.0)]]),
        ]
        mp = GI.MultiPolygon(polys)
        env = Extents.Extent(X = (0.0, 10.0), Y = (0.0, 10.0))
        sss = GO.extract_segment_strings(relate_geom(mp), GO.GEOM_A, env)
        @test length(sss) == 1
        @test sss[1].pts[1] == (0.0, 0.0)
        # ring-level filtering: shell intersects the extent but the hole does not
        poly_far_hole = GI.Polygon([
            [(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (0.0, 100.0), (0.0, 0.0)],
            [(90.0, 90.0), (90.0, 95.0), (95.0, 95.0), (95.0, 90.0), (90.0, 90.0)],
        ])
        sss = GO.extract_segment_strings(relate_geom(poly_far_hole), GO.GEOM_A, env)
        @test length(sss) == 1
        @test sss[1].ring_id == 0
    end

    @testset "line extraction removes repeated points" begin
        line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (2.0, 2.0)])
        sss = GO.extract_segment_strings(relate_geom(line), GO.GEOM_A, nothing)
        @test length(sss) == 1
        ss = only(sss)
        @test ss.dim == GO.DIM_L
        @test ss.ring_id == -1
        @test ss.parent_polygonal === nothing
        @test ss.pts == [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)]
        @test !GO.is_closed(ss)
    end

    # Regression (found by cross-validation, Task 24): geometries backed by
    # StaticArrays — e.g. `GO.extent_to_polygon` output, whose ring wraps an
    # SVector — must extract to plain `Vector` point lists; a typed
    # comprehension over `GI.getpoint` collected them to a `SizedVector`,
    # which `_orient_ring` rejected.
    @testset "static-array-backed ring extraction" begin
        poly = GO.extent_to_polygon(Extents.Extent(X = (0.0, 1.0), Y = (0.0, 1.0)))
        sss = GO.extract_segment_strings(relate_geom(poly), GO.GEOM_A, nothing)
        @test length(sss) == 1
        @test only(sss).pts isa Vector{Tuple{Float64, Float64}}
        @test GO.relate_predicate(GO.RelateNG(), GO.pred_equalstopo(), poly, poly)
    end

    # CW shell ring: (0 0, 0 10, 10 10, 10 0, 0 0) after reorientation
    shell = only(GO.extract_segment_strings(
        relate_geom(GI.Polygon([
            [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
        ])), GO.GEOM_A, nothing))
    open_line = only(GO.extract_segment_strings(
        relate_geom(GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])),
        GO.GEOM_A, nothing))

    @testset "prev/next vertex with ring wraparound" begin
        @test shell.pts == [(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)]
        # node at ring closure vertex: prev wraps to last distinct vertex
        @test GO.prev_vertex(shell, 1, (0.0, 0.0)) == (10.0, 0.0)
        @test GO.next_vertex(shell, 4, (0.0, 0.0)) == (0.0, 10.0)
        # node interior to a segment: prev/next are the segment endpoints
        @test GO.prev_vertex(shell, 1, (0.0, 5.0)) == (0.0, 0.0)
        @test GO.next_vertex(shell, 1, (0.0, 5.0)) == (0.0, 10.0)
        # node at a non-closure vertex
        @test GO.prev_vertex(shell, 2, (0.0, 10.0)) == (0.0, 0.0)
        @test GO.next_vertex(shell, 1, (0.0, 10.0)) == (10.0, 10.0)
        # open line endpoints have no prev/next beyond the ends
        @test GO.prev_vertex(open_line, 1, (0.0, 0.0)) === nothing
        @test GO.next_vertex(open_line, 2, (2.0, 2.0)) === nothing
    end

    @testset "create_node_section" begin
        # vertex node at the ring closure point
        ns = GO.create_node_section(shell, 1, GO.vertex_node((0.0, 0.0)))
        @test ns.is_a
        @test ns.dim == GO.DIM_A
        @test ns.id == 1 && ns.ring_id == 0
        @test ns.is_node_at_vertex
        @test ns.node == GO.vertex_node((0.0, 0.0))
        @test GO.get_vertex(ns, 0) == (10.0, 0.0)
        @test GO.get_vertex(ns, 1) == (0.0, 10.0)
        # proper crossing node: section vertices are the segment endpoints
        node = GO.crossing_node((0.0, 0.0), (0.0, 10.0), (-1.0, 5.0), (1.0, 5.0))
        ns = GO.create_node_section(shell, 1, node)
        @test !ns.is_node_at_vertex
        @test ns.node == node
        @test GO.get_vertex(ns, 0) == (0.0, 0.0)
        @test GO.get_vertex(ns, 1) == (0.0, 10.0)
        # vertex node at the end of an open line: no next vertex
        ns = GO.create_node_section(open_line, 2, GO.vertex_node((2.0, 2.0)))
        @test ns.is_node_at_vertex
        @test ns.dim == GO.DIM_L && ns.ring_id == -1
        @test GO.get_vertex(ns, 0) == (1.0, 1.0)
        @test GO.get_vertex(ns, 1) === nothing
    end

    @testset "is_containing_segment" begin
        # at segment start vertex - always contained
        @test GO.is_containing_segment(shell, 1, (0.0, 0.0))
        # in segment interior - always contained
        @test GO.is_containing_segment(shell, 1, (0.0, 5.0))
        # at segment end vertex - assigned to the next segment
        @test !GO.is_containing_segment(shell, 1, (0.0, 10.0))
        # closed ring: final segment endpoint belongs to the first segment
        @test !GO.is_containing_segment(shell, 4, (0.0, 0.0))
        # open line: final segment contains its endpoint
        @test GO.is_containing_segment(open_line, 2, (2.0, 2.0))
        @test !GO.is_containing_segment(open_line, 1, (1.0, 1.0))
    end
end
