using Test
 
import GeoInterface as GI, GeometryOps as GO
using GeoInterface: Extents

poly = GI.Polygon([GI.LinearRing([(1, 2), (3, 4), (5, 6), (1, 2)]), 
                   GI.LinearRing([(3, 4), (5, 6), (6, 7), (3, 4)])])

@test_all_implementations "embed_extent" poly begin
    ext_poly = GO.embed_extent(poly)
    lr1, lr2 = GI.getgeom(ext_poly)
    @test ext_poly.extent == Extents.Extent(X=(1, 6), Y=(2, 7))
    @test lr1.extent == Extents.Extent(X=(1, 5), Y=(2, 6))
    @test lr2.extent == Extents.Extent(X=(3, 6), Y=(4, 7))
end
