using Test

import GeoInterface as GI
using GeometryOpsTestHelpers

@testset "JTS XML fixture reader" begin
    xml = """
    <run>
      <desc>Small synthetic relate/overlay fixture</desc>
      <precisionModel scale="1.0" offsetx="2.5" offsety="-3.5"/>
      <case>
        <desc>PP equal</desc>
        <a>POINT (0 0)</a>
        <b>POINT (0 0)</b>
        <test><op name="relate" arg1="A" arg2="B" arg3="0FFFFFFF2">true</op></test>
        <test><op name="intersection" arg1="A" arg2="B">POINT (0 0)</op></test>
      </case>
    </run>
    """

    mktemp() do path, io
        write(io, xml)
        close(io)

        test_set = load_test_set(path)
        @test length(test_set.cases) == 1
        @test test_set.description == "Small synthetic relate/overlay fixture"
        @test test_set.precision_model isa JTSPrecisionModel
        @test test_set.precision_model.scale == 1.0
        @test test_set.precision_model.offsetx == 2.5
        @test test_set.precision_model.offsety == -3.5
        @test fixture_family("TestRelatePP.xml") == :relate
        @test fixture_family("TestNGOverlayP.xml") == :overlay

        case = only(test_set.cases)
        @test case.description == "PP equal"
        @test geometry_category(case.geom_a) == :point
        @test case_category(case) == :point_point
        @test primary_conformance_category(test_set, case) == :point_point
        @test conformance_categories(test_set, case) == [:point_point, :precision_snap]
        @test has_conformance_category(test_set, case, :point_point)
        @test length(case.operations) == 2

        relate = case.operations[1]
        @test is_relate_operation(relate)
        @test is_runnable(relate)
        @test relate.argument_refs == ["A", "B", "0FFFFFFF2"]
        @test relate.arguments[1] === case.geom_a
        @test relate.arguments[2] === case.geom_b
        @test relate.arguments[3] == "0FFFFFFF2"
        @test relate.expected === true

        overlay = case.operations[2]
        @test is_overlay_operation(overlay)
        @test is_overlay_operation("intersectionNG")
        @test is_overlay_operation("symDifferenceNG")
        @test is_relate_operation("coveredBy")
        @test is_relate_operation("equalsTopo")
        @test GI.trait(overlay.expected) isa GI.PointTrait

        overlay_cases = load_test_cases(path; operations = ["intersection"])
        @test length(overlay_cases) == 1
        @test length(only(overlay_cases).operations) == 1
        @test only(only(overlay_cases).operations).name == "intersection"

        marked = load_test_set(
            path;
            broken = JTSFixtureRule(operation = "relate", reason = "matrix later"),
            unimplemented = JTSFixtureRule(
                case = r"PP equal",
                operation = r"intersection",
                reason = "overlay later",
            ),
        )
        marked_ops = only(marked.cases).operations
        @test is_broken(marked_ops[1])
        @test marked_ops[1].status_reason == "matrix later"
        @test is_unimplemented(marked_ops[2])
        @test marked_ops[2].status_reason == "overlay later"
    end
end

@testset "JTS conformance categories" begin
    mktempdir() do dir
        relate_dir = joinpath(dir, "robust")
        mkpath(relate_dir)
        relate_path = joinpath(relate_dir, "TestRobustRelateSynthetic.xml")
        overlay_path = joinpath(dir, "TestNGOverlayGCSynthetic.xml")

        write(
            relate_path,
            """
            <run>
              <precisionModel type="FIXED" scale="10.0"/>
              <case>
                <desc>zero length line</desc>
                <a>LINESTRING (0 0, 0 0)</a>
                <b>LINESTRING (0 0, 1 1)</b>
                <test><op name="intersects" arg1="A" arg2="B">true</op></test>
              </case>
              <case>
                <desc>polygon with hole</desc>
                <a>POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 2 1, 2 2, 1 1))</a>
                <b>POLYGON ((3 3, 5 3, 5 5, 3 5, 3 3))</b>
                <test><op name="relate" arg1="A" arg2="B" arg3="*********">true</op></test>
              </case>
              <case>
                <desc>relate collection</desc>
                <a>GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 1))</a>
                <b>POINT (0 0)</b>
                <test><op name="intersects" arg1="A" arg2="B">true</op></test>
              </case>
            </run>
            """,
        )

        write(
            overlay_path,
            """
            <run>
              <case>
                <desc>overlay collection</desc>
                <a>GEOMETRYCOLLECTION (POINT (0 0), POINT (1 1))</a>
                <b>POINT (0 0)</b>
                <test><op name="intersection" arg1="A" arg2="B">POINT (0 0)</op></test>
              </case>
              <case>
                <desc>mixed overlay collection</desc>
                <a>GEOMETRYCOLLECTION (POINT (0 0), LINESTRING (0 0, 1 1))</a>
                <b>POINT (0 0)</b>
                <test><op name="intersection" arg1="A" arg2="B">POINT (0 0)</op></test>
              </case>
            </run>
            """,
        )

        relate_set = load_test_set(relate_path)
        zero_line, hole_case, relate_collection = relate_set.cases

        @test conformance_categories(relate_set, zero_line) == [
            :line_line,
            :zero_length_line,
            :repeated_coordinates,
            :precision_snap,
            :robust_failure,
        ]
        @test conformance_categories(relate_set, hole_case) == [
            :area_area,
            :holes_and_touching_rings,
            :precision_snap,
            :robust_failure,
        ]
        @test primary_conformance_category(relate_set, relate_collection) ==
              :relateng_collection

        overlay_set = load_test_set(overlay_path)
        overlay_collection, mixed_overlay_collection = overlay_set.cases
        @test primary_conformance_category(overlay_set, overlay_collection) ==
              :overlayng_collection
        @test primary_conformance_category(overlay_set, mixed_overlay_collection) == :other

        zero_length_only = load_test_set(relate_path; categories = :zero_length_line)
        @test length(zero_length_only.cases) == 1
        @test only(zero_length_only.cases).description == "zero length line"

        inventory = conformance_inventory([relate_set, overlay_set])
        @test inventory[:line_line] == 1
        @test inventory[:area_area] == 1
        @test inventory[:relateng_collection] == 1
        @test inventory[:overlayng_collection] == 1
        @test inventory[:other] == 1
        @test inventory[:precision_snap] == 3
        @test inventory[:robust_failure] == 3
    end
end

@testset "JTS empty WKT placeholders" begin
    empty_point = jts_wkt_to_geom("POINT EMPTY")
    @test empty_point isa JTSEmptyGeometry
    @test geometry_category(empty_point) == :empty

    raw_collection = jts_wkt_to_geom(
        "GEOMETRYCOLLECTION (POLYGON EMPTY, LINESTRING (0 0, 1 1))",
    )
    @test raw_collection isa JTSRawGeometry
    @test geometry_category(raw_collection) == :collection

    raw_wkb = jts_wkt_to_geom("0101000000000000000000F03F0000000000000040")
    @test raw_wkb isa JTSRawGeometry
    @test geometry_category(raw_wkb) == :unknown

    raw_linearring = jts_wkt_to_geom("LINEARRING (0 0, 1 0, 1 1, 0 0)")
    @test raw_linearring isa JTSRawGeometry
    @test geometry_category(raw_linearring) == :line
end

@testset "JTS XML fixture discovery and batch loading" begin
    mktempdir() do dir
        relate_path = joinpath(dir, "TestRelateSynthetic.xml")
        overlay_path = joinpath(dir, "TestNGOverlaySynthetic.xml")
        write(
            relate_path,
            """
            <run>
              <case>
                <desc>raw input is still loadable</desc>
                <a>010100000000000000000000000000000000000000</a>
                <b>POINT (0 0)</b>
                <test><op name="intersects" arg1="A" arg2="B">false</op></test>
              </case>
            </run>
            """,
        )
        write(
            overlay_path,
            """
            <run>
              <case>
                <a>POINT (0 0)</a>
                <b>POINT (1 1)</b>
                <test><op name="union" arg1="A" arg2="B">MULTIPOINT ((0 0), (1 1))</op></test>
              </case>
            </run>
            """,
        )

        @test find_jts_test_files(dir) == sort([overlay_path, relate_path])
        @test find_jts_test_files(dir; family = :relate) == [relate_path]
        @test find_jts_test_files(dir; filename = r"Overlay") == [overlay_path]

        test_sets = load_test_sets(dir; operations = r"intersects")
        @test length(test_sets) == 2
        relate_set = only(filter(ts -> basename(ts.filepath) == "TestRelateSynthetic.xml", test_sets))
        overlay_set = only(filter(ts -> basename(ts.filepath) == "TestNGOverlaySynthetic.xml", test_sets))
        @test only(relate_set.cases).geom_a isa JTSRawGeometry
        @test only(relate_set.cases).operations[1].expected === false
        @test isempty(overlay_set.cases)

        skipped = load_test_set(
            relate_path;
            skip = JTSFixtureRule(file = "TestRelateSynthetic", operation = :intersects),
        )
        @test is_skipped(only(only(skipped.cases).operations))
    end
end
