using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS
using GeometryOpsTestHelpers

# All triangle vertices lie on the square boundary, but the closing edge leaves it.
# Shared vertices do not establish containment.
@testset "Shared-boundary enclave" begin
    notch_points = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (2.0, 2.0)]
    fill_points = [(0.0, 0.0), (2.0, 2.0), (0.0, 4.0)]
    alg = GO.FosterHormannClipping(GO.Planar())
    cache = GO.FosterHormannCache(alg)
    closed_polygon(points) = GI.Polygon([vcat(points, [first(points)])])
    for shift in 0:2, reverse_notch in (false, true), reverse_fill in (false, true)
        np = circshift(notch_points, shift)
        fp = circshift(fill_points, shift)
        notch = closed_polygon(reverse_notch ? reverse(np) : np)
        fill = closed_polygon(reverse_fill ? reverse(fp) : fp)
        for (a, b) in ((notch, fill), (fill, notch))
            @test_implementations isempty(GO.intersection(alg, $a, $b; target=GI.PolygonTrait()))
            @test_implementations GO.intersection_area(alg, $a, $b) == 0.0
            @test GO.intersection_area(alg, a, b; cache) == 0.0
            @test_implementations sum(GO.area, GO.union(alg, $a, $b; target=GI.PolygonTrait())) == 16.0
            @test_implementations sum(GO.area, GO.difference(alg, $a, $b; target=GI.PolygonTrait())) == GO.area($a)
        end
        # An entirely shared boundary still denotes coincident polygons.
        @test_implementations GO.intersection_area(alg, $notch, $notch) == 12.0
    end
    # All triangle vertices also lie on a square, this time with edges inside.
    square = closed_polygon([(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
    triangle = closed_polygon([(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)])
    for (a,b) in ((square,triangle), (triangle,square))
        @test_implementations GO.intersection_area(alg, $a, $b) == 8.0
        @test_implementations sum(GO.area, GO.intersection(alg, $a, $b; target=GI.PolygonTrait())) == 8.0
    end
end
