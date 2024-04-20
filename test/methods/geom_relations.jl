using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

# Tests of DE-9IM Methods
pt1 = LG.Point([0.0, 0.0])
pt2 = LG.Point([5.0, 5.0])
pt3 = LG.Point([1.0, 0.0])
pt4 = LG.Point([0.5, 0.0])
pt5 = LG.Point([0.5, 0.25])
pt6 = LG.Point([0.6, 0.4])
pt7 = LG.Point([0.4, 0.8])

l1 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
l2 = LG.LineString([[0.0, 0.0], [1.0, 0.0]])
l3 = LG.LineString([[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
l4 = LG.LineString([[0.5, 0.0], [1.0, 0.0], [1.0, 0.5]])
l5 = LG.LineString([[0.0, 0.0], [-1.0, -1.0]])
l6 = LG.LineString([[2.0, 2.0], [0.0, 1.0]])
l7 = LG.LineString([[0.5, 1.0], [0.5, -1.0]])
l8 = LG.LineString([[0.0, 0.0], [0.5, 0.0], [0.5, 0.5], [1.0, -0.5]])
l9 = LG.LineString([[0.0, 1.0], [0.0, -1.0], [1.0, 1.0]])
l10 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]])
l11 = LG.LineString([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0], [-1.0, 0.0]])
l12 = LG.LineString([[0.6, 0.5], [0.6, 0.9]])
l13 = LG.LineString([[2.0, 2.0], [3.0, 3.0]])
l14 = LG.LineString([[0.6, 0.25], [0.6, 0.35]])
l15 = LG.LineString([[0.0, 3.0], [4.0, 3.0]])
l16 = LG.LineString([[0.3, -0.7], [1.0, 0.0], [3.0, 0.6]])

r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]])
r2 = LG.LinearRing([[0.5, 0.2], [0.6, 0.4], [0.7, 0.2], [0.5, 0.2]])
r3 = LG.LinearRing([[0.2, 0.7], [0.4, 0.9], [0.5, 0.6], [0.2, 0.7]])
r4 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]])
r5 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 1.0], [1.0, 2.0], [1.0, 1.0], [0.0, 0.0]])
r6 = LG.LinearRing([[0.0, 0.0], [-1.0, 0.0], [-1.0, 1.0], [0.0, 1.0], [0.0, 0.0]])
r7 = LG.LinearRing([[0.5, 0.5], [1.5, 0.5], [1.5, 1.5], [0.5, 1.5], [0.5, 0.5]])
r8 = LG.LinearRing([[0.1, 0.1], [0.2, 0.1], [0.2, 0.2], [0.1, 0.2], [0.1, 0.1]])
r9 = LG.LinearRing([[1.0, -0.5], [0.5, -1.0], [1.5, -1.0], [1.0, -0.5]])
r10 = LG.LinearRing([[0.5, 0.2], [0.6, 0.4], [0.3, 0.3], [0.5, 0.2]])
r11 = LG.LinearRing([[0.55, 0.21], [0.55, 0.23], [0.65, 0.23], [0.66, 0.21], [0.55, 0.21]])

p1 = LG.Polygon(r4, [r2, r3])
p2 = LG.Polygon(r4)
p3 = LG.Polygon(r2)
p4 = LG.Polygon(r6)
p5 = LG.Polygon(r9)
p6 = LG.Polygon(r11)
p7 = LG.Polygon(r7)
p8 = LG.Polygon([[[0.0, 0.0], [0.0, 3.0], [1.0, 2.0], [2.0, 3.0], [3.0, 2.0], [4.0, 3.0], [4.0, 0.0], [0.0, 0.0]]])
p9 = LG.Polygon([[[0.0, 0.0], [0.0, 3.0], [1.0, 4.0], [2.0, 3.0], [3.0, 4.0], [4.0, 3.0], [4.0, 0.0], [0.0, 0.0]]])
p10 = LG.Polygon([
    [[0.1, 0.5], [0.1, 0.99], [0.6, 0.99], [0.6, 0.5], [0.1, 0.5]],
    [[0.15, 0.55], [0.15, 0.95], [0.55, 0.95], [0.55, 0.55], [0.15, 0.55]]
])
p11 = LG.Polygon(r3)

mpt1 = LG.MultiPoint([pt1, pt2])
mpt2 = LG.MultiPoint([pt2, pt3])
mpt3 = LG.MultiPoint([pt4, pt5])
ml1 = LG.MultiLineString([l5, l6, l7])
ml2 = LG.MultiLineString([l1])
mp1 = LG.MultiPolygon([p1])
mp2 = LG.MultiPolygon([p6, p7])
gc1 = LG.GeometryCollection([pt1, l5, p6])

test_pairs = [
    # Points and geometries
    (pt1, pt1, "pt1", "pt1", "Same point"),
    (pt1, pt2, "pt1", "pt2", "Different point"),
    (pt1, l1, "pt1", "l1", "Point on line endpoint"),
    (pt2, l1, "pt2", "l1", "Point outside line"),
    (pt3, l1, "pt3", "l1", "Point on line segment"),
    (pt4, l1, "pt4", "l1", "Point on line vertex between segments"),
    (l1, pt3, "l1", "pt3", "Point on line segment (order swapped)"),
    (pt1, r1, "pt1", "r1", "Point on ring 'endpoint'"),
    (pt2, r1, "pt2", "r1", "Point outside ring"),
    (pt3, r1, "pt3", "r1", "Point on ring segment"),
    (pt4, r1, "pt4", "r1", "Point on ring vertex between segments"),
    (r1, pt3, "r1", "pt3", "Point on ring segment (order swapped)"),
    (p1, pt1, "p1", "pt1", "Point on vertex of polygon"),
    (pt2, p1, "pt2", "p1", "Point outside of polygon's external ring"),
    (pt4, p1, "pt4", "p1", "Point on polygon's edge"),
    (pt5, p1, "pt5", "p1", "Point inside of polygon"),
    (pt6, p1, "pt6", "p1", "Point on hole edge"),
    (pt7, p1, "pt7", "p1", "Point inside of polygon hole"),
    (p1, pt5, "p1", "pt5", "Point inside of polygon (order swapped)"),
    # # Lines and geometries
    (l1, l1, "l1", "l1", "Same line"),
    (l2, l1, "l2", "l1", "L2 is one segment of l1"),
    (l3, l1, "l3", "l1", "L3 shares one segment with l1 and has one segment outside"),
    (l4, l1, "l4", "l1", "L4 shares half of each of l1 segments"),
    (l5, l1, "l5", "l1", "L5 shares one endpoint with l1 but not segments"),
    (l6, l1, "l6", "l1", "Lines are disjoint"),
    (l7, l1, "l7", "l1", "L7 crosses through one of l1's segments"),
    (l8, l1, "l8", "l1", "Overlaps one segment and crosses through another segment"),
    (l9, l1, "l9", "l1", "Two segments touch and two segments crosses"),
    (l16, l1, "l16", "l1", "L16 bounces off of l1's corner"),
    (l1, r1, "l1", "r1", "Line inside of ring"),
    (l3, r1, "l3", "r1", "Line covers one edge of linear ring and has segment outside"),
    (l5, r1, "l5", "r1", "Line and linear ring are only covered by vertex"),
    (l6, r1, "l6", "r1", "Line and ring are disjoint"),
    (l7, r1, "l7", "r1", "Line crosses through two ring edges"),
    (l8, r1, "l8", "r1", "Line crosses through two ring edges and touches third edge"),
    (l10, r1, "l10", "r1", "Line is equal to linear ring"),
    (l11, r1, "l11", "r1", "Line covers linear ring and then has extra segment"),
    (l1, p1, "l1", "p1", "Line on polygon edge"),
    (l3, p1, "l3", "p1", "Line on polygon edge and extending beyond polygon edge"),
    (l5, p1, "l5", "p1", "Line outside polygon connected by a vertex"),
    (l7, p1, "l7", "p1", "Line through polygon cutting to the outside"),
    (l12, p1, "l12", "p1", "Line inside polygon"),
    (l13, p1, "l13", "p1", "Line outside of polygon"),
    (l14, p1, "l14", "p1", "Line in polygon hole"),
    (l15, p8, "l15", "p8", "Line outside crown-shaped polygon but touching edges"),
    (l15, p9, "l15", "p9", "Line within crown-shaped polygon but touching edges"),
    # Ring and geometries
    (r1, l1, "r1", "l1", "Line is within linear ring"),
    (r1, l3, "r1", "l3", "Line covers one edge of linear ring and has segment outside"),
    (r1, l5, "r1", "l5", "Line and linear ring are only connected at vertex"),
    (r1, l6, "r1", "l6", "Line and linear ring are disjoint"),
    (r1, l7, "r1", "l7", "Line crosses though two ring edges"),
    (r1, l8, "r1", "l8", "Line crosses through two ring edges and touches third edge"),
    (r1, l10, "r1", "l10", "Line is equal to linear ring"),
    (r1, l11, "r1", "l11", "Line covers linear ring and then has extra segment"),
    (r1, r1, "r1", "r1", "Same rings"), 
    (r2, r1, "r2", "r1", "Disjoint ring with one 'inside' of hole created"),
    (r3, r1, "r3", "r1", "Disjoint ring with one 'outside' of hole created"),
    (r4, r1, "r4", "r1", "Rings share two sides and rest of sides dont touch"),
    (r1, r5, "r1", "r5", "Ring shares all edges with other ring, plus an extra loop"),
    (r1, r6, "r1", "r6", "Rings share just one vertex"),
    (r1, r7, "r1", "r7", "Rings cross one another"),
    (r4, p1, "r4", "p1", "Ring on boundary of polygon"),
    (r1, p1, "r1", "p1", "Ring on boundary and cutting through polygon"),
    (r2, p1, "r2", "p1", "Ring on hole bounday"),
    (r6, p1, "r6", "p1", "Ring touches polygon at one vertex"),
    (r7, p1, "r7", "p1", "Ring crosses through polygon"),
    (r8, p1, "r8", "p1", "Ring inside of polygon"),
    (r9, p1, "r9", "p1", "Ring outside of polygon"),
    (r10, p1, "r10", "p1", "Ring inside of polygon and shares hole's edge"),
    (r11, p1, "r11", "p1", "Ring inside of polygon hole"),
    # Polygon and geometries
    (p1, p1, "p1", "p1", "Same polygons"),
    (p1, p2, "p1", "p2", "P1 and p2 are the same but p1 has holes"),
    (p2, p1, "p2", "p1", "P1 and p2 are the same but p1 has holes (order swapped)"),
    (p3, p1, "p3", "p1", "P3 is equal to one of p1's holes"),
    (p4, p1, "p4", "p1", "Polygon's share just one vertex"),
    (p5, p1, "p5", "p1", "Polygon outside of other polygon"),
    (p6, p1, "p6", "p1", "Polygon inside of other polygon's hole"),
    (p7, p1, "p7", "p1", "Polygons overlap"),
    (p10, p1, "p10", "p1", "Polygon's with nested holes"),
    # Multigeometries
    (mpt1, mpt1, "mpt1", "mpt1", "Same set of points for multipoints"),
    (mpt1, mpt2, "mpt1", "mpt2", "Some point matches, others are different"),
    (mpt1, mpt3, "mpt1", "mpt3", "No shared points"),
    (ml1, ml2, "ml1", "ml2", "Lines in ml1 cross and touch ml2"),
    (mp1, mp2, "mp1", "mp2", "Polygons in mp1 are inside hole and overlap"),
    (gc1, ml1, "gc1", "ml1", "Make sure collection works with multi-geom"),
]

function test_geom_relation(GO_f, LG_f, f_name; swap_points = false)
    for (g1, g2, sg1, sg2, sdesc) in test_pairs
        if swap_points
            g1, g2 = g2, g1
            sg1, sg2 = sg2, sg1
        end
        go_val = GO_f(g1, g2)
        lg_val = LG_f(g1, g2)
        @test go_val == lg_val
        go_val != lg_val && println("\n↑ TEST INFO: $sg1 $f_name $sg2 - $sdesc \n\n")
    end
end

@testset "Contains" begin test_geom_relation(GO.contains, LG.contains, "contains"; swap_points = true) end
@testset "Covered By" begin test_geom_relation(GO.coveredby, LG.coveredby, "coveredby") end
@testset "Covers" begin test_geom_relation(GO.covers, LG.covers, "covers"; swap_points = true) end
@testset "Disjoint" begin test_geom_relation(GO.disjoint, LG.disjoint, "disjoint")end
@testset "Intersect" begin test_geom_relation(GO.intersects, LG.intersects, "intersects") end
@testset "Touches" begin test_geom_relation(GO.touches, LG.touches, "touches") end
@testset "Within" begin test_geom_relation(GO.within, LG.within, "within") end


@testset "Overlaps" begin
    @testset "Points/MultiPoints" begin
        p1 = LG.Point([0.0, 0.0])
        p2 = LG.Point([0.0, 1.0])
        # Two points can't overlap
        @test GO.overlaps(p1, p1) == LG.overlaps(p1, p2)
    
        mp1 = LG.MultiPoint([[0.0, 1.0], [4.0, 4.0]])
        mp2 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0]])
        mp3 = LG.MultiPoint([[0.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
        # No shared points, doesn't overlap
        @test GO.overlaps(p1, mp1) == LG.overlaps(p1, mp1)
        # One shared point, does overlap
        @test GO.overlaps(p2, mp1) == LG.overlaps(p2, mp1)
        # All shared points, doesn't overlap
        @test GO.overlaps(mp1, mp1) == LG.overlaps(mp1, mp1)
        # Not all shared points, overlaps
        @test GO.overlaps(mp1, mp2) == LG.overlaps(mp1, mp2)
        # One set of points entirely inside other set, doesn't overlap
        @test GO.overlaps(mp2, mp3) == LG.overlaps(mp2, mp3)
        # Not all points shared, overlaps
        @test GO.overlaps(mp1, mp3) == LG.overlaps(mp1, mp3)
    
        mp1 = LG.MultiPoint([
            [-36.05712890625, 26.480407161007275],
            [-35.7220458984375, 27.137368359795584],
            [-35.13427734375, 26.83387451505858],
            [-35.4638671875, 27.254629577800063],
            [-35.5462646484375, 26.86328062676624],
            [-35.3924560546875, 26.504988828743404],
        ])
        mp2 = GI.MultiPoint([
            [-35.4638671875, 27.254629577800063],
            [-35.5462646484375, 26.86328062676624],
            [-35.3924560546875, 26.504988828743404],
            [-35.2001953125, 26.12091815959972],
            [-34.9969482421875, 26.455820238459893],
        ])
        # Some shared points, overlaps
        @test GO.overlaps(mp1, mp2) == LG.overlaps(mp1, mp2)
        @test GO.overlaps(mp1, mp2) == GO.overlaps(mp2, mp1)
    end
    
    @testset "Lines/Rings" begin
        l1 = LG.LineString([[0.0, 0.0], [0.0, 10.0]])
        l2 = LG.LineString([[0.0, -10.0], [0.0, 20.0]])
        l3 = LG.LineString([[0.0, -10.0], [0.0, 3.0]])
        l4 = LG.LineString([[5.0, -5.0], [5.0, 5.0]])
        # Line can't overlap with itself
        @test GO.overlaps(l1, l1) == LG.overlaps(l1, l1)
        # Line completely within other line doesn't overlap
        @test GO.overlaps(l1, l2) == GO.overlaps(l2, l1) == LG.overlaps(l1, l2)
        # Overlapping lines
        @test GO.overlaps(l1, l3) == GO.overlaps(l3, l1) == LG.overlaps(l1, l3)
        # Lines that don't touch
        @test GO.overlaps(l1, l4) == LG.overlaps(l1, l4)
        # Linear rings that intersect but don't overlap
        r1 = LG.LinearRing([[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]])
        r2 = LG.LinearRing([[1.0, 1.0], [1.0, 6.0], [6.0, 6.0], [6.0, 1.0], [1.0, 1.0]])
        @test GO.overlaps(r1, r2) == LG.overlaps(r1, r2)
    end
    
    @testset "Polygons/MultiPolygons" begin
        p1 = LG.Polygon([[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]])
        p2 = LG.Polygon([
            [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
            [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
        ])
        # Test basic polygons that don't overlap
        @test GO.overlaps(p1, p2) == LG.overlaps(p1, p2)
        @test !GO.overlaps(p1, (1, 1))
        @test !GO.overlaps((1, 1), p2)
    
        p3 = LG.Polygon([[[1.0, 1.0], [1.0, 6.0], [6.0, 6.0], [6.0, 1.0], [1.0, 1.0]]])
        # Test basic polygons that overlap
        @test GO.overlaps(p1, p3) == LG.overlaps(p1, p3)
    
        p4 = LG.Polygon([[[20.0, 5.0], [20.0, 10.0], [18.0, 10.0], [18.0, 5.0], [20.0, 5.0]]])
        # Test one polygon within the other
        @test GO.overlaps(p2, p4) == GO.overlaps(p4, p2) == LG.overlaps(p2, p4)
    
        p5 = LG.Polygon(
            [[
                [-53.57208251953125, 28.287451910503744],
                [-53.33038330078125, 28.29228897739706],
                [-53.34136352890625, 28.430052892335723],
                [-53.57208251953125, 28.287451910503744],
            ]]
        )
        # Test equal polygons
        @test GO.overlaps(p5, p5) == LG.overlaps(p5, p5)
    
        # Test multipolygons
        m1 = LG.MultiPolygon([
            [[[0.0, 0.0], [0.0, 5.0], [5.0, 5.0], [5.0, 0.0], [0.0, 0.0]]],
            [
                [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
                [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]]
            ]
        ])
        # Test polygon that overlaps with multipolygon
        @test GO.overlaps(m1, p3) == LG.overlaps(m1, p3)
        # Test polygon in hole of multipolygon, doesn't overlap
        @test GO.overlaps(m1, p4) == LG.overlaps(m1, p4)
    end
end
@testset "Crosses" begin
	line6 = GI.LineString([(1.0, 1.0), (1.0, 2.0), (1.0, 3.0), (1.0, 4.0)])
	poly7 = GI.Polygon([[(-1.0, 2.0), (3.0, 2.0), (3.0, 3.0), (-1.0, 3.0), (-1.0, 2.0)]])

	@test GO.crosses(GI.LineString([(-2.0, 2.0), (4.0, 2.0)]), line6) == true
	@test GO.crosses(GI.LineString([(0.5, 2.5), (1.0, 1.0)]), poly7) == true
	@test GO.crosses(GI.MultiPoint([(1.0, 2.0), (12.0, 12.0)]), GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])) == true
	@test GO.crosses(GI.MultiPoint([(1.0, 0.0), (12.0, 12.0)]), GI.LineString([(1, 1), (1, 2), (1, 3), (1, 4)])) == false
	@test GO.crosses(GI.LineString([(-2.0, 2.0), (-4.0, 2.0)]), poly7) == false
end