using JSON3
import GeoFormatTypes as GFT, GeoInterface as GI, GeometryOps as GO
import TGGeometry
import LibGEOS
import WellKnownGeometry


using Test

TG_PRED_SYMBOL_TO_FUNCTION = Dict(
    [sym => getproperty(GO, sym) for sym in TGGeometry.TG_PREDICATES]
)

TG_IGNORE_LIST = Set([:crosses, :overlaps, :equals])

testsets = JSON3.read(read(joinpath(@__DIR__, "data", "mpoly.jsonc")))

# # assume GEOS is the ultimate source of truth

testset = testsets[4]
geoms = GO.tuples.(GFT.WellKnownText.((GFT.Geom(),), testset.geoms))



function run_testsets(json_file_path)
    testsets = JSON3.read(read(json_file_path))
    for i in eachindex(testsets)
        testset = testsets[i]
        geoms = GO.tuples.(GFT.WellKnownText.((GFT.Geom(),), testset.geoms))
        @testset let testset_index = i
            for (predname, results) in testset.predicates
                !haskey(TG_PRED_SYMBOL_TO_FUNCTION, predname) && continue
                predname in TG_IGNORE_LIST && continue
                predicate_f = TG_PRED_SYMBOL_TO_FUNCTION[predname]

                @testset let predicate = predname
                    go_result = try
                        predicate_f(geoms[1], geoms[2])
                    catch e
                        println("Error: $e")
                        @test_broken false
                        continue
                    end
                    lg_result = predicate_f(GO.GEOS(), geoms[1], geoms[2])
                    # tg_result = predicate_f(GO.TG(), geoms[1], geoms[2])

                    expected = first(results) == "T"

                    @test go_result == expected
                    @test lg_result == expected
                    # @test tg_result == expected
                end
            end
        end
    end
end

@testset "All TG testsets" begin
    for file in filter(endswith(".jsonc"), readdir(joinpath(@__DIR__, "data"); join = true))
        filename = splitext(basename(file))[1]
        @testset "$filename" begin
            run_testsets(file)
        end
    end
end
