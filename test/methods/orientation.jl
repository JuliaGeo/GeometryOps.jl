using Test, GeometryOps
import GeoInterface as GI
import GeometryOps as GO

line1 = GI.LineString([[9.170356, 45.477985], [9.164434, 45.482551], [9.166644, 45.484003]])
line2 = GI.LineString([[9.169356, 45.477985], [9.163434, 45.482551], [9.165644, 45.484003]])
line3 = GI.LineString([
    (-111.544189453125, 24.186847428521244),
    (-110.687255859375, 24.966140159912975),
    (-110.4510498046875, 24.467150664739002),
    (-109.9951171875, 25.180087808990645)
])
line4 = GI.LineString([
    (-111.4617919921875, 24.05148034322011),
    (-110.8795166015625, 24.681961205014595),
    (-110.841064453125, 24.14174098050432),
    (-109.97863769531249, 24.617057340809524)
])

poly1 = GI.Polygon([[[0, 0], [1, 0], [1, 1], [0.5, 0.5], [0, 1], [0, 0]]])
poly2 = GI.Polygon([[[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]]])

l1 = GI.LineString([[0, 0], [1, 1], [1, 0], [0, 0]])
l2 = GI.LineString([[0, 0], [1, 0], [1, 1], [0, 0]])

@test_all_implementations "Orientation" (poly1, poly2, l1, l2) begin

    # @test isparallel(line1, line2) == true
    # @test isparallel(line3, line4) == false
    
    @test isconcave(poly1) == true
    @test isconcave(poly2) == false

    @test isclockwise(l1) == true
    @test isclockwise(l2) == false
end
