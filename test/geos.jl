import GeometryOps as GO, GeoInterface as GI, LibGEOS as LG
using GeometryOps: GEOS

pt1 = GI.Point((0.0, 0.0))
mpt1 = GI.MultiPoint([pt1, pt1])
l1 = GI.Line([(0.0, 0.0), (0.0, 1.0)])

concave_coords = [(0.0, 0.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]
l2 = GI.LineString(concave_coords)
l3 = GI.LineString(concave_coords[1:(end - 1)])
r1 = GI.LinearRing(concave_coords)
r2 = GI.LinearRing(concave_coords[1:(end - 1)])
r3 = GI.LinearRing([(1.0, 1.0), (1.0, 1.5), (1.5, 1.5), (1.5, 1.0), (1.0, 1.0)])
concave_angles = [90.0, 270.0, 90.0, 90.0, 90.0, 90.0]

p1 = GI.Polygon([r3])
p2 = GI.Polygon([[(0.0, 0.0), (0.0, 4.0), (3.0, 0.0), (0.0, 0.0)]])
p3 = GI.Polygon([[(-3.0, -2.0), (0.0,0.0), (5.0, 0.0), (-3.0, -2.0)]])
p4 = GI.Polygon([r1])
p5 = GI.Polygon([r1, r3])

mp1 = GI.MultiPolygon([p2, p3])
c1 = GI.GeometryCollection([pt1, l2, p2])

@testset "GeometryOpsLibGEOSExt with GeoInterface Geometries" begin
    @testset "Functionality Tests" begin
        @testset "Buffer" begin
            @test GO.buffer(p1, 1.0) == LG.buffer(p1, 1.0)
        end
        # @testset "DE-9IM"
        #     function test_geom_relation(GO_f, LG_f, f_name; swap_points = false)
        #         for (g1, g2, sg1, sg2, sdesc) in test_pairs
        #             if swap_points
        #                 g1, g2 = g2, g1
        #                 sg1, sg2 = sg2, sg1
        #             end
        #             go_val = GO_f(g1, g2)
        #             lg_val = LG_f(g1, g2)
        #             @test go_val == lg_val
        #             go_val != lg_val && println("\n↑ TEST INFO: $sg1 $f_name $sg2 - $sdesc \n\n")
        #         end
        #     end
            
        #     @testset "Contains" begin test_geom_relation(GO.contains, LG.contains, "contains"; swap_points = true) end
        #     @testset "Covered By" begin test_geom_relation(GO.coveredby, LG.coveredby, "coveredby") end
        #     @testset "Covers" begin test_geom_relation(GO.covers, LG.covers, "covers"; swap_points = true) end
        #     @testset "Disjoint" begin test_geom_relation(GO.disjoint, LG.disjoint, "disjoint")end
        #     @testset "Intersect" begin test_geom_relation(GO.intersects, LG.intersects, "intersects") end
        #     @testset "Touches" begin test_geom_relation(GO.touches, LG.touches, "touches") end
        #     @testset "Within" begin test_geom_relation(GO.within, LG.within, "within") end
        # end
    end
end
