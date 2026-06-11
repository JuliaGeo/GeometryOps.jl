# Executable overlay-test loop, moved verbatim from `jts_testset_reader.jl`
# (which is now a pure parser). Not wired into CI; run manually.
include(joinpath(@__DIR__, "jts_testset_reader.jl"))

# TODO: parameterize; path predates vendored test/data/jts files
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
