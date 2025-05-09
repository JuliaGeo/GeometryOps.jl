using Test
using DataFrames
import GeometryOps as GO, GeoInterface as GI

using NaturalEarth

@testset "DataFrames extension: can we use copycols=false?" begin
    df = DataFrame(naturalearth("admin_0_countries", 110); copycols = true)
    transformed = GO.transform(identity, df; copycols = true)
    transformed_lazy = GO.transform(identity, df; copycols = false)

    @test transformed.NAME == transformed_lazy.NAME
    @test transformed.NAME !== df.NAME # test that they are different arrays
    @test transformed_lazy.NAME === df.NAME # test that they are the same array
end