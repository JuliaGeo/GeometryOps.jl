using Test

import GeoInterface as GI
using GeometryOpsTestHelpers

@testset "JTS XML fixture reader" begin
    xml = """
    <run>
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
        @test fixture_family("TestRelatePP.xml") == :relate
        @test fixture_family("TestNGOverlayP.xml") == :overlay

        case = only(test_set.cases)
        @test case.description == "PP equal"
        @test geometry_category(case.geom_a) == :point
        @test case_category(case) == :point_point
        @test length(case.operations) == 2

        relate = case.operations[1]
        @test is_relate_operation(relate)
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
end
