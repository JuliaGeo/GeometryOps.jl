@testset "ArchGDAL" begin
    # Yoinked from ArchGDAL tests
    p1 = AG.fromWKT("POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))")
    @test GI.area(p1) ≈ abs(GO.signed_area(p1))
    @test GO.signed_area(p1) > 0 # test that the signed area is positive
end