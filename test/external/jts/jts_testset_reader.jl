using XML

import WellKnownGeometry
import GeoFormatTypes as GFT
import GeometryOps as GO
import GeoInterface as GI

"""
    jts_wkt_to_geom(wkt::String)

Convert a JTS WKT string to a GeometryOps geometry, via WellKnownGeometry.jl and GO.tuples.

The reason this exists is because WellKnownGeometry doesn't work well with newlines in subsidiary geometries,
so this sanitizes the input before parsing and converting.
"""
function jts_wkt_to_geom(wkt::String)
    sanitized_wkt = join(strip.(split(wkt, "\n")), "")
    geom = GFT.WellKnownText(GFT.Geom(), sanitized_wkt)
    return GO.tuples(geom)
end

struct TestItem{T}
    operation::String
    arg1::GO.GI.Wrappers.WrapperGeometry
    arg2::GO.GI.Wrappers.WrapperGeometry
    expected_result::T
end

Base.show(io::IO, ::MIME"text/plain", item::TestItem) = print(io, "TestItem(operation = $(item.operation), expects $(GI.trait(item.expected_result)))")
Base.show(io::IO, item::TestItem) = show(io, MIME"text/plain"(), item)

struct Case
    description::String
    geom_a::GO.GI.Wrappers.WrapperGeometry
    geom_b::GO.GI.Wrappers.WrapperGeometry
    items::Vector{TestItem}
end

function load_test_cases(filepath::String)
    doc = read(filepath, XML.Node) # lazy parsing
    run = only(children(doc))
    test_cases = Case[]
    for case in children(run)
        if tag(case) != "case"
            continue
        end
        push!(test_cases, parse_case(case))
    end
    return test_cases
end

function parse_case(case::XML.Node)
    description = value(only(children(case.children[1])))
    a = jts_wkt_to_geom(value(only(children(case.children[2]))))
    b = jts_wkt_to_geom(value(only(children(case.children[3]))))

    items = TestItem[]
    for item in children(case)[4:end]
        ops = children(item)
        for op in ops
            op_attrs = XML.attributes(op)
            operation = op_attrs["name"]
            arg1 = op_attrs["arg1"]
            arg2 = op_attrs["arg2"]
            expected_result = jts_wkt_to_geom(value(only(op.children)))
            push!(items, TestItem(operation, lowercase(arg1) == "a" ? a : b, lowercase(arg2) == "a" ? a : b, expected_result))
        end
    end
    return Case(description, a, b, items)
end


testfile = "/Users/anshul/.julia/dev/geo/jts/modules/tests/src/test/resources/testxml/general/TestOverlayAA.xml"

cases = load_test_cases(testfile)

using Test

for case in cases
    @testset "$(case.description)" begin
        for item in case.items
            @testset "$(item.operation)" begin
                result = if item.operation == "union"
                     GO.union(item.arg1, item.arg2; target = GO.PolygonTrait())
                elseif item.operation == "difference"
                    GO.difference(item.arg1, item.arg2; target = GO.PolygonTrait())
                elseif item.operation == "intersection"
                    GO.intersection(item.arg1, item.arg2; target = GO.PolygonTrait())
                elseif item.operation == "symdifference"
                    continue
                else
                    continue
                end

                finalresult = if length(result) == 0
                    nothing
                elseif length(result) == 1
                    only(result)
                else
                    GI.MultiPolygon(result)
                end

                if isnothing(finalresult)
                    @warn("No result")
                    continue
                end

                if GI.geomtrait(item.expected_result) isa Union{GI.MultiPolygonTrait, GI.PolygonTrait}
                    difference_in_areas = GO.area(GO.difference(finalresult, item.expected_result; target = GO.PolygonTrait()))
                    if difference_in_areas > 1
                        @warn("Difference in areas: $(difference_in_areas)")
                        f, a, p = poly(finalresult; label = "Final result (GO)", axis = (; title = "$(case.description) - $(item.operation)"))
                        poly!(a, item.expected_result; label = "Expected result (JTS)")
                        display(f)
                    end
                    # @test difference_in_areas < 1
                else
                    @test_broken GO.equals(finalresult, item.expected_result)
                end
            end
        end
    end
end