using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
using ..TestHelpers

# Test of polygon clipping
p1 = GI.Polygon([[(0.0, 0.0), (5.0, 5.0), (10.0, 0.0), (5.0, -5.0), (0.0, 0.0)]])
p2 = GI.Polygon([[(3.0, 0.0), (8.0, 5.0), (13.0, 0.0), (8.0, -5.0), (3.0, 0.0)]])
p3 = GI.Polygon([[(4.526700198111509, 3.4853728532584696), (2.630732683726619, -4.126134282323841),
    (-0.7638522032421201, -4.418734350277446), (-4.367920073785058, -0.2962672719707883),
    (4.526700198111509, 3.4853728532584696)]])
p4 = GI.Polygon([[(5.895141140952208, -0.8095078714426418), (2.8634927670695283, -4.625511746720306),
    (-1.790623183259246, -4.138092164660989), (-3.9856656502985843, -0.5275687876429914),
    (-2.554809853598822, 3.553455552936806), (1.1865909598835922, 4.984203644564732),
    (5.895141140952208, -0.8095078714426418)]])
p5 = GI.Polygon([[(2.5404227081738795, 0.5995497066446837), (0.7215593458353178, -1.5811392990170074),
    (-0.6792714151561866, -1.0909218208298457), (-0.5721092724334685, 2.0387826195734795),
    (0.0011462224659918308, 2.3273077404755487), (2.5404227081738795, 0.5995497066446837)]])
p6 = GI.Polygon([[(3.2022522653586183, -4.4613815131276615), (-1.0482425878695998, -4.579816661708281),
    (-3.630239248625253, 2.0443933767558677), (-2.6164940041615927, 3.4708149011067224),
    (1.725945294696213, 4.954192017601067), (3.2022522653586183, -4.4613815131276615)]])
p7 = GI.Polygon([[(1.2938349167338743, -3.175128530227131), (-2.073885870841754, -1.6247711001754137),
    (-5.787437985975053, 0.06570713422599561), (-2.1308128111898093, 5.426689675486368),
    (2.3058074184797244, 6.926652158268195), (1.2938349167338743, -3.175128530227131)]])
p8 = GI.Polygon([[(-2.1902469793743924, -1.9576242117579579), (-4.726006206053999, 1.3907098941556428),
    (-3.165301985923147, 2.847612825874245), (-2.5529280962099428, 4.395492123980911),
    (0.5677700216973937, 6.344638314896882), (3.982554842356183, 4.853519613487035),
    (5.251193948893394, 0.9343031382106848), (5.53045582244555, -3.0101433691361734),
    (-2.1902469793743924, -1.9576242117579579)]])
p9 = GI.Polygon([[(0.0, 0.0), (0.0, 4.0), (7.0, 4.0), (7.0, 0.0), (0.0, 0.0)]])
p10 = GI.Polygon([[(1.0, -3.0), (1.0, 1.0), (3.5, -1.5), (6.0, 1.0), (6.0, -3.0), (1.0, -3.0)]])
p11 = GI.Polygon([[(1.0, 1.0), (4.0, 1.0), (4.0, 2.0), (2.0, 2.0), (2.0, 3.0), (4.0, 3.0),
    (4.0, 4.0), (1.0, 4.0), (1.0, 1.0)]])
p12 = GI.Polygon([[(3.0, 0.0), (5.0, 0.0), (5.0, 5.0), (3.0, 5.0), (3.0, 0.0)]])
p13 = GI.Polygon([[(1.0, 1.0), (4.0, 1.0), (4.0, 2.0), (2.0, 2.0), (2.0, 3.0), (4.0, 3.0),
    (4.0, 4.0), (2.0, 4.0), (2.0, 5.0), (4.0, 5.0), (4.0, 6.0), (2.0, 6.0), (2.0, 7.0),
    (4.0, 7.0), (4.0, 8.0), (1.0, 8.0), (1.0, 1.0)]])
p14 = GI.Polygon([[(3.0, 0.0), (5.0, 0.0), (5.0, 9.0), (3.0, 9.0), (3.0, 0.0)]])
p15 = GI.Polygon([[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]])
p16 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p17 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0), (0.0, 0.0)], [(2.0, 1.0), (3.0, 1.0), (3.0, 2.0), (2.0, 2.0), (2.0, 1.0)]])
p18 = GI.Polygon([[(1.0, -1.0), (1.0, 4.0), (5.0, 4.0), (5.0, -1.0), (1.0, -1.0)]])
p19 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)], [(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])
p20 = GI.Polygon([[(2.0, -1.0), (5.0, -1.0), (5.0, 5.0), (2.0, 5.0), (2.0, -1.0)]])
p21 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p22 = GI.Polygon([[(1.0, -1.0), (2.0, -1.0), (2.0, 4.0), (1.0, 4.0), (1.0, -1.0)]])
p23 = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 7.0), (0.0, 7.0), (0.0, 0.0)],
    [(1.0, 1.0), (5.0, 1.0), (5.0, 6.0), (1.0, 6.0), (1.0, 1.0)]])
p24 = GI.Polygon([[(2.0, 2.0), (8.0, 2.0), (8.0, 5.0), (2.0, 5.0), (2.0, 2.0)],
    [(3.0, 3.0), (7.0, 3.0), (7.0, 4.0), (3.0, 4.0), (3.0, 3.0)]])
p25 = GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 8.0), (0.0, 8.0), (0.0, 0.0)],
    [(4.0, 0.5), (4.5, 0.5), (4.5, 3.5), (4.0, 3.5), (4.0, 0.5)],
    [(2.0, 4.0), (4.0, 4.0), (4.0, 6.0), (2.0, 6.0), (2.0, 4.0)]])
p26 = GI.Polygon([[(3.0, 1.0), (8.0, 1.0), (8.0, 7.0), (3.0, 7.0), (3.0, 5.0), (6.0, 5.0), (6.0, 3.0), (3.0, 3.0), (3.0, 1.0)],
    [(3.5, 5.5), (6.0, 5.5), (6.0, 6.5), (3.5, 6.5), (3.5, 5.5)],
    [(5.5, 1.5), (5.5, 2.5), (3.5, 2.5), (3.5, 1.5), (5.5, 1.5)]])
p27 = GI.Polygon([[(0.0, 0.0), (8.0, 0.0), (10.0, -1.0), (8.0, 1.0), (8.0, 2.0), (7.0, 5.0), (6.0, 4.0), (3.0, 5.0), (3.0, 3.0), (0.0, 0.0)]])
p28 = GI.Polygon([[(1.0, 1.0), (3.0, -1.0), (6.0, 2.0), (8.0, 0.0), (8.0, 4.0), (4.0, 4.0), (1.0, 1.0)]])
p29 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (3.0, 1.0), (1.0, 1.0), (0.0, 2.0), (0.0, 0.0)]])
p30 = GI.Polygon([[(4.0, 0.0), (3.0, 1.0), (1.0, 1.0), (0.0, 0.0), (0.0, 2.0), (4.0, 2.0), (4.0, 0.0)]])
p31 = GI.Polygon([[(0.0, 0.0), (2.0, 1.0), (4.0, 0.0), (2.0, 4.0), (1.0, 2.0), (0.0, 3.0), (0.0, 0.0)]])
p32 = GI.Polygon([[(4.0, 3.0), (3.0, 2.0), (4.0, 2.0), (4.0, 3.0)]])
p33 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p34 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 0.0)]])
p35 = GI.Polygon([[(1.0, 0.0), (2.0, 0.0), (2.0, -1.0), (1.0, -1.0), (1.0, 0.0)]])
p36 = GI.Polygon([[(2.0, 1.0), (3.0, 0.0), (4.0, 1.0), (3.0, 3.0), (2.0, 1.0)]])
p37 = GI.Polygon([[(1.0, -1.0), (2.0, -1.0), (2.0, -2.0), (1.0, -2.0), (1.0, -1.0)]])
p38 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)],
    [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]])
p39 = GI.Polygon([[(5.0, 0.0), (8.0, 0.0), (8.0, 3.0), (5.0, 3.0), (5.0, 0.0)],
    [(6.0, 1.0), (7.0, 1.0), (7.0, 2.0), (7.0, 1.0), (6.0, 1.0)]])
p40 = GI.Polygon([[(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0), (0.0, 0.0)],
    [(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)]])
p41 = GI.Polygon([[(3.0, -1.0), (10.0, -1.0), (10.0, 9.0), (3.0, 9.0), (3.0, -1.0)],
    [(4.0, 3.0), (5.0, 3.0), (5.0, 4.0), (4.0, 4.0), (4.0, 3.0)],
    [(6.0, 1.0), (7.0, 1.0), (7.0, 2.0), (6.0, 2.0), (6.0, 1.0)]])
p42 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
    [(1.0, 1.0), (3.0, 1.0), (3.0, 1.5), (1.0, 1.5), (1.0, 1.0)],
    [(1.0, 2.5), (3.0, 2.5), (3.0, 3.0), (1.0, 3.0), (1.0, 2.5)]])
p43 = GI.Polygon([[(2.0, -1.0), (5.0, -1.0), (5.0, 5.0), (2.0, 5.0), (2.0, -1.0)],
    [(3.5, 4.0), (4.5, 4.0), (4.5, 4.5), (3.5, 4.5), (3.5, 4.0)],
    [(3.5, 3.0), (4.5, 3.0), (4.5, 3.5), (3.5, 3.5), (3.5, 3.0)],
    [(3.5, 2.0), (4.5, 2.0), (4.5, 2.5), (3.5, 2.5), (3.5, 2.0)]])
p44 = GI.Polygon([[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)],
    [(5.0, 2.0), (5.0, 7.0), (7.0, 7.0), (7.0, 2.0), (5.0, 2.0)],
    [(8.0, 2.0), (8.0, 7.0), (9.0, 7.0), (9.0, 2.0), (8.0, 2.0)],
    [(7.25, 3.5), (7.25, 3.75), (7.75, 3.75), (7.75, 3.5), (7.25, 3.5)]])
p45 = GI.Polygon([[(3.0, -2.0), (3.0, 8.0), (13.0, 8.0), (13.0, -2.0), (3.0, -2.0)],
    [(5.5, 2.5), (5.5, 3.0), (8.5, 3.0), (8.5, 2.5), (5.5, 2.5)],
    [(5.5, 4.0), (5.5, 4.5), (8.5, 4.5), (8.5, 4.0), (5.5, 4.0)]])
p46 = GI.Polygon([[(0.0, 0.0), (6.0, 0.0), (6.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p47 = GI.Polygon([[(1.0, -1.0), (1.0, 1.0), (2.0, 0.0), (3.0, 0.0), (3.0, -1.0), (1.0, -1.0)]])
p48 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (3.0, 0.0), (5.0, 0.0), (5.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p49 = GI.Polygon([[(1.0, -1.0), (4.0, -1.0), (4.0, 2.0), (3.0, 2.0), (3.0, 0.0), (2.0, 1.0), (1.0, 0.0), (1.0, -1.0)]])
p50 = GI.Polygon([[(0.0, 0.0), (1.5, 0.0), (1.5, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 3.0), (0.0, 0.0)]])
p51 = GI.Polygon([[(1.0, 1.0), (2.0, 1.0), (2.0, -1.0), (1.5, -1.0), (1.5, -1.0), (1.0, -1.0), (1.0, 1.0), (1.0, 1.0)]])
p52 = GI.Polygon([[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)],
    [(1.0, 8.0), (2.0, 8.0), (2.0, 9.0), (1.0, 9.0), (1.0, 8.0)],
    [(2.5, 2.5), (2.5, 7.5), (7.5, 7.5), (7.5, 2.5), (2.5, 2.5)],
    [(8.0, 1.0), (8.0, 2.0), (9.0, 2.0), (9.0, 1.0), (8.0, 1.0)]])
p53 = GI.Polygon([[(3.0, 3.0), (3.0, 6.0), (6.0, 6.0), (6.0, 3.0), (3.0, 3.0)],
    [(4.0, 4.0), (4.0, 5.0), (5.0, 5.0), (5.0, 4.0), (4.0, 4.0)]])
p54 = GI.Polygon([[(2.5, 2.5), (2.5, 7.5), (7.5, 7.5), (7.5, 2.5), (2.5, 2.5)],
    [(4.0, 4.0), (4.0, 5.0), (5.0, 5.0), (5.0, 4.0), (4.0, 4.0)]])
p55 = GI.Polygon([[(5.0, 0.25), (5.0, 5.0), (9.5, 5.0), (9.5, 0.25), (5.0, 0.25)],
    [(6.0, 3.0), (6.0, 4.0), (7.0, 4.0), (7.0, 3.0), (6.0, 3.0)],
    [(7.5, 0.5), (7.5, 2.5), (9.25, 2.5), (9.25, 0.5), (7.5, 0.5)]])

mp1 = GI.MultiPolygon([p1])
mp2 = GI.MultiPolygon([p2])
mp3 = GI.MultiPolygon([p34, p36])
mp4 = GI.MultiPolygon([p33, p35])
    
test_pairs = [
    (p1, p1, "p1", "p1", "Same polygon"),
    (p1, p2, "p1", "p2", "Convex polygons that intersect (diamonds, four vertices)"),
    (p3, p4, "p3", "p4", "Convex polygons that intersect (randomly generated, many edges)"),
    (p5, p6, "p5", "p6", "Convex polygons that intersect (randomly generated, many edges)"),
    (p7, p8, "p7", "p8", "Concave polygons that intersect (randomly generated, many edges)"),
    (p5, p7, "p5", "p7", "Convex and concave polygons intersect"),
    (p9, p10, "p9", "p10", "Figure 10 from Greiner Hormann paper"),
    (p11, p12, "p11", "p12", "Polygons whose intersection is two distinct regions"),
    (p13, p14, "p13", "p14", "Polygons whose intersection is two distinct regions"),
    (p15, p16, "p15", "p16", "First polygon in second polygon"),
    (p16, p15, "p16", "p15", "Second polygon in first polygon"),
    (p17, p18, "p17", "p18", "First polygon with a hole (hole completely in second polygon), second without a hole"),
    (p18, p17, "p18", "p17", "First polygon with no hole, second with a hole (hole completely in first polygon)"),
    (p19, p20, "p19", "p20", "First polygon with a hole (hole not completely in second polygon), second without a hole"),
    (p42, p20, "p42", "p20", "First polygon with two holes (holes not completely in second polygon), second without a hole"),
    (p20, p19, "p20", "p19", "First polygon with no hole, second with a hole (hole not completely in first polygon)"),
    (p20, p42, "p20", "p42", "First polygon with no holes, second with two holes (holes not completely in second polygon)"),
    (p21, p22, "p21", "p22", "Polygons form cross, splitting each other"),
    (p23, p24, "p23", "p24", "Polygons are both donuts with intersecting holes"),
    (p25, p26, "p25", "p26", "Polygons both have two holes that intersect in various ways"),
    (p27, p28, "p27", "p28", "Figure 12 from Foster extension for degeneracies"),
    (p29, p30, "p29", "p30", "Figure 13 from Foster extension for degeneracies"),
    (p31, p32, "p31", "p32", "Polygons touch at just one point"),
    (p33, p34, "p33", "p34", "One polygon inside of the other, sharing an edge"),
    (p33, p35, "p33", "p35", "Polygons outside of one another, sharing an edge"),
    (p33, p36, "p33", "p36", "All intersection points are V-intersections as defined by Foster"),
    (p33, p37, "p33", "p37", "Polygons are completely disjoint (no holes)"),
    (p38, p39, "p38", "p39", "Polygons are completely disjoint (both have one hole)"),
    (p40, p41, "p40", "p41", "Two overlapping polygons with three total holes in overlap region"),
    (p42, p43, "p42", "p43", "First polygon 2 holes, second polygon 3 holes. Holes do not overlap"),
    (p43, p42, "p43", "p42", "First polygon 3 holes, second polygon 2 holes. Holes do not overlap"),
    (p44, p45, "p44", "p45", "Holes form a ring, with an additional hole within that ring of holes"),
    (p46, p47, "p46", "p47", "Intersecting polygons that share one edge and cross through the same edge"),
    (p48, p49, "p48", "p49", "Intersecting polygons that share two edges in a row and cross through the different edge"),
    (p50, p51, "p50", "p51", "Intersection polygons with opposite winding orders and repeated points"),
    (p52, p53, "p52", "p53", "Polygon with two holes completely inside of one (of three) holes of other polygon"),
    (p52, p54, "p52", "p54", "Polygon with two holes has exterior equal to one (of three) holes of other polygon"),
    (p52, p55, "p52", "p55", "Polygon within another polygon, with intersecting and disjoint holes"),
    (p1, mp2, "p1", "mp2", "Polygon overlapping with mutlipolygon with one sub-polygon"),
    (mp1, p2, "mp1", "p2", "Polygon overlapping with mutlipolygon with one sub-polygon"),
    (mp1, mp2, "mp1", "mp2", "Two overlapping multipolygons, which each have one sub-polygon"),
    (mp3, p33, "mp3", "p33", "Mutipolygon, with two sub-polygon, both of which overlap with other polygon"),
    (p33, mp3, "p33", "mp3", "Polygon overlaps with both sub-polygons of a multipolygon"),
    (mp3, p34, "mp3", "p34", "Multipolygon, with two sub-polygon, overlaps with polygon equal to sub-polygon"),
    (p34, mp3, "p34", "mp3", "Polygon overlaps with multipolygon, where on of the sub-polygons is equivalent"),
    (mp3, p35, "mp3", "p35", "Mulitpolygon where just one sub-polygon touches other polygon"),
    (p35, mp3, "p35", "mp3", "Polygon that touches just one of the sub-polygons of multipolygon"),
    (mp4, mp3, "mp4", "mp3", "Two multipolygons, which with two sub-polygons, where some of the sub-polygons intersect and some don't")
]

const ϵ = 1e-10
# Compare clipping results from GeometryOps and LibGEOS
function compare_GO_LG_clipping(GO_f, LG_f, p1, p2)
    
    LG_result_geom = LG_f(p1, p2)
    if LG_result_geom isa LG.GeometryCollection
        poly_list = LG.Polygon[]
        for g in GI.getgeom(LG_result_geom)
            g isa LG.Polygon && push!(poly_list, g)
        end
        LG_result_geom = LG.MultiPolygon(poly_list)
    end

    for _accelerator in (GO.AutoAccelerator(), GO.NestedLoop(), GO.SingleSTRtree(), GO.SingleNaturalTree(), #=GO.DoubleNaturalTree(), =# #=GO.ThinnedDoubleNaturalTree(), =# #=GO.DoubleSTRtree()=#)
    @testset let accelerator = _accelerator # this is a ContextTestSet that is otherwise invisible but adds context to the testset
        GO_result_list = GO_f(GO.FosterHormannClipping(accelerator), p1, p2; target = GI.PolygonTrait())
        # Check if nothing is returned
        if isempty(GO_result_list) && (LG.isEmpty(LG_result_geom) || LG.area(LG_result_geom) == 0)
            @test true
            continue
        end
        # Check for unnecessary points
        @test !(sum(GI.npoint, GO_result_list; init = 0.0) > GI.npoint(LG_result_geom))
        # Make sure last point is repeated
        for poly in GO_result_list
            for ring in GI.getring(poly)
                @test !(GI.getpoint(ring, 1) != GI.getpoint(ring, GI.npoint(ring)))
            end
        end

        # Check if polygons cover the same area
        local GO_result_geom
        if length(GO_result_list) == 1
            GO_result_geom = GO_result_list[1]
        else
            GO_result_geom = GI.MultiPolygon(GO_result_list)
        end
        
        diff_1_poly = @test_nowarn LG.difference(GO_result_geom, LG_result_geom)
        diff_2_poly = @test_nowarn LG.difference(LG_result_geom, GO_result_geom)
        if GI.isempty(diff_1_poly) || GI.isempty(diff_2_poly)
            @test true # if both differences are empty, then the areas are the same by default.
        else
            diff_1_area = LG.area(diff_1_poly)
            diff_2_area = LG.area(diff_2_poly)
            @test diff_1_area ≤ ϵ && diff_2_area ≤ ϵ
        end
    end # testset
    end # loop
end

# Test clipping functions and print error message if tests fail
function test_clipping(GO_f, LG_f, f_name)
    for (p1, p2, sg1, sg2, sdesc) in test_pairs
        @testset_implementations "$sg1 $f_name $sg2 - $sdesc" begin
            compare_GO_LG_clipping(GO_f, LG_f, $p1, $p2) # this executes tests internally
        end
    end
    return
end

@testset "Intersection" begin test_clipping(GO.intersection, LG.intersection, "intersection") end
@testset "Union" begin test_clipping(GO.union, LG.union, "union") end
@testset "Difference" begin test_clipping(GO.difference, LG.difference, "difference") end

@testset "TracingError show" begin
    message = "Test tracing error"
    poly_a_small = p1
    poly_b_small = p2

    poly_a_large, poly_b_large = GO.segmentize.((poly_a_small, poly_b_small), max_distance = 0.001)

    # PolyNode now has two type parameters: {T, P} where P is the point type
    PT = Tuple{Float64, Float64}
    empty_nodes = GO.PolyNode{Float64, PT}[]

    # Test that the message is contained in the exception
    @test_throws "Test tracing error" throw(GO.TracingError(message, poly_a_small, poly_b_small, empty_nodes, empty_nodes, Int[]))
    # Test that the coordinates of the polygons are contained in the exception,
    # if the polygon is small enough
    @test_throws "$(GI.coordinates(poly_a_small))" throw(GO.TracingError(message, poly_a_small, poly_b_small, empty_nodes, empty_nodes, Int[]))
    # Test that the coordinates of the polygons are not contained in the exception,
    # if the polygon is large enoughs
    @test_throws "The polygons are contained in the exception object" throw(GO.TracingError(message, poly_a_large, poly_b_large, empty_nodes, empty_nodes, Int[]))
end

@testset "Spherical polygon clipping pipeline" begin
    # Test that the spherical clipping pipeline runs without error
    # Note: Full spherical clipping requires additional work on the tracing algorithm
    # to properly handle spherical geometry. Currently, intersection point detection
    # works correctly but polygon tracing still uses planar assumptions.

    # Two overlapping squares on the sphere (small enough to be well-behaved)
    poly1 = GI.Polygon([[
        (-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0), (-10.0, -10.0)
    ]])
    poly2 = GI.Polygon([[
        (0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0), (0.0, 0.0)
    ]])

    # Verify the pipeline runs without throwing
    result = @test_nowarn GO.intersection(GO.Spherical(), poly1, poly2; target=GI.PolygonTrait())
    @test result isa Vector

    # Verify intersection_points works correctly with Spherical manifold
    # This is the part that's fully implemented
    line1 = GI.LineString([(-10.0, 0.0), (10.0, 0.0)])
    line2 = GI.LineString([(0.0, -10.0), (0.0, 10.0)])
    pts = GO.intersection_points(GO.Spherical(), line1, line2)
    @test length(pts) == 1
    @test pts[1] isa GO.UnitSpherical.UnitSphericalPoint
    # Convert to geographic to check longitude/latitude
    to_geo = GO.UnitSpherical.GeographicFromUnitSphere()
    geo = to_geo(pts[1])
    @test isapprox(geo[1], 0.0, atol=1e-6)  # longitude
    @test isapprox(geo[2], 0.0, atol=1e-6)  # latitude

    # Test that spherical area function works on results
    if !isempty(result)
        area_result = @test_nowarn GO.area(GO.Spherical(), result[1])
        @test area_result isa Number
    end
end

@testset "Spherical intersection vs union correctness" begin
    # This test verifies that intersection and union are not swapped.
    # Intersection should produce smaller area, union should produce larger area.

    # Two overlapping polygons on the sphere
    poly1 = GI.Polygon([[
        (-15.0, -15.0), (15.0, -15.0), (15.0, 15.0), (-15.0, 15.0), (-15.0, -15.0)
    ]])
    poly2 = GI.Polygon([[
        (0.0, 0.0), (30.0, 0.0), (30.0, 30.0), (0.0, 30.0), (0.0, 0.0)
    ]])

    # Compute areas of input polygons
    area1 = GO.area(GO.Spherical(), poly1)
    area2 = GO.area(GO.Spherical(), poly2)
    min_input_area = min(area1, area2)
    max_input_area = max(area1, area2)

    # Compute intersection
    intersection_result = GO.intersection(GO.Spherical(), poly1, poly2; target=GI.PolygonTrait())
    @test !isempty(intersection_result)
    intersection_area = sum(p -> GO.area(GO.Spherical(), p), intersection_result)

    # Compute union
    union_result = GO.union(GO.Spherical(), poly1, poly2; target=GI.PolygonTrait())
    @test !isempty(union_result)
    union_area = sum(p -> GO.area(GO.Spherical(), p), union_result)

    # Key invariants that would fail if intersection/union are swapped:
    # 1. Intersection area must be smaller than the smallest input polygon
    @test intersection_area < min_input_area
    # 2. Union area must be larger than the largest input polygon
    @test union_area > max_input_area
    # 3. Intersection area must be smaller than union area
    @test intersection_area < union_area
    # 4. Inclusion-exclusion: area(A) + area(B) = area(A ∪ B) + area(A ∩ B)
    @test isapprox(area1 + area2, union_area + intersection_area, rtol=0.01)
end

@testset "Spherical intersection - vertex on edge cases" begin
    # Test case where a vertex of polygon A lies ON an edge of polygon B
    # (not at a vertex of B). This previously caused a BoundsError due to
    # duplicate intersection nodes being created.

    # Case 1: Vertex of A on meridian edge of B
    # poly_a has a vertex at (70, -80) which lies on poly_b's edge from
    # (70, -79.5) to (70, -80.5) - a meridian line at longitude 70
    poly_a = GI.Polygon([[(70.0, -80.0), (70.0, -79.0), (71.0, -79.0), (71.0, -80.0), (70.0, -80.0)]])
    poly_b = GI.Polygon([[(69.0, -80.5), (69.0, -79.5), (70.0, -79.5), (70.0, -80.5), (69.0, -80.5)]])

    # This should not throw - previously threw BoundsError accessing index 0
    result1 = @test_nowarn GO.intersection(GO.Spherical(), poly_a, poly_b; target=GI.PolygonTrait())
    # The polygons share only a vertex touch, so intersection should be empty or a point
    @test length(result1) == 0  # No area overlap, just vertex touch

    # Case 2: Two rectangles with shared meridian edge
    # poly_c and poly_d share an edge along longitude 70
    poly_c = GI.Polygon([[(70.0, -80.0), (70.0, -79.0), (71.0, -79.0), (71.0, -80.0), (70.0, -80.0)]])
    poly_d = GI.Polygon([[(69.0, -80.0), (69.0, -79.0), (70.0, -79.0), (70.0, -80.0), (69.0, -80.0)]])

    # Should handle shared edge correctly
    result2 = @test_nowarn GO.intersection(GO.Spherical(), poly_c, poly_d; target=GI.PolygonTrait())
    # Shared edge but no area overlap
    @test length(result2) == 0

    # Case 3: Overlapping rectangles (should produce non-empty intersection)
    poly_e = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    poly_f = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

    result3 = @test_nowarn GO.intersection(GO.Spherical(), poly_e, poly_f; target=GI.PolygonTrait())
    @test length(result3) == 1  # Should have one intersection polygon

    # Case 4: Polygons near poles (high latitude stress test)
    # This tests numerical precision in spherical calculations
    poly_g = GI.Polygon([[(160.0, 38.0), (160.0, 39.0), (161.0, 39.0), (161.0, 38.0), (160.0, 38.0)]])
    poly_h = GI.Polygon([[(160.0, 38.24), (160.0, 39.19), (161.2, 39.19), (161.2, 38.24), (160.0, 38.24)]])

    result4 = @test_nowarn GO.intersection(GO.Spherical(), poly_g, poly_h; target=GI.PolygonTrait())
    @test length(result4) == 1  # Should have one intersection polygon
end

@testset "Spherical intersection - UnitSphericalPoint input" begin
    # Test that intersection works correctly when polygons are already
    # in UnitSphericalPoint coordinates (not lat/long).
    # Previously this failed because _tuple_point extracted only (x, y) from
    # 3D Cartesian coordinates, and UnitSphereFromGeographic then wrongly
    # interpreted these as lat/long.

    # Create overlapping lat/long polygons
    poly_a_ll = GI.Polygon([[(69.0, -81.0), (69.0, -80.0), (70.0, -80.0), (70.0, -81.0), (69.0, -81.0)]])
    poly_b_ll = GI.Polygon([[(69.0, -80.5), (69.0, -79.5), (70.0, -79.5), (70.0, -80.5), (69.0, -80.5)]])

    # Convert to UnitSphericalPoint
    transform = GO.UnitSpherical.UnitSphereFromGeographic()
    poly_a_sph = GO.transform(transform, poly_a_ll)
    poly_b_sph = GO.transform(transform, poly_b_ll)

    # Both paths should give the same result
    result_ll = GO.intersection(GO.Spherical(), poly_a_ll, poly_b_ll; target=GI.PolygonTrait())
    result_sph = GO.intersection(GO.Spherical(), poly_a_sph, poly_b_sph; target=GI.PolygonTrait())

    @test length(result_ll) == length(result_sph)
    @test length(result_ll) == 1  # These polygons overlap

    # Check that the intersection areas are approximately equal
    area_ll = GO.area(GO.Spherical(), result_ll[1])
    area_sph = GO.area(GO.Spherical(), result_sph[1])
    @test isapprox(area_ll, area_sph, rtol=0.01)
end

@testset "Spherical intersection - vertex on interior of edge (α≈1 case)" begin
    # Test case where a vertex of polygon A lies ON an edge of polygon B,
    # but NOT at an endpoint of that edge. This creates an intersection at
    # α ≈ 1 (end of an A edge) that must be properly tracked and applied
    # to the next vertex.
    #
    # Previously, intersections at α ≈ 1 were skipped, expecting them to be
    # caught at α ≈ 0 of the next edge. But spherical_arc_intersection only
    # returns the interior crossing, not the vertex-on-edge intersection.

    using GeometryOps.UnitSpherical: UnitSphericalPoint

    # Two overlapping quadrilaterals where vertices of A lie on edges of B
    # (and vice versa), but not at endpoints
    p1_pts = [
        UnitSphericalPoint(0.3019858622586644, 0.8296993375502143, -0.46947156278589075),
        UnitSphericalPoint(0.3047421791070115, 0.8372722558519103, -0.4539904997395468),
        UnitSphericalPoint(0.29008334968929284, 0.8424632196002839, -0.4539904997395468),
        UnitSphericalPoint(0.2874596182894698, 0.8348433503288182, -0.46947156278589075),
        UnitSphericalPoint(0.3019858622586644, 0.8296993375502143, -0.46947156278589075)
    ]

    p2_pts = [
        UnitSphericalPoint(0.30139703746681035, 0.8280815547305803, -0.47269521313549934),
        UnitSphericalPoint(0.303920062408707, 0.8350135087871621, -0.4586774856131918),
        UnitSphericalPoint(0.2895744973775592, 0.840098592125015, -0.45867304910807816),
        UnitSphericalPoint(0.28716095284187626, 0.8331277227671596, -0.47269100342586656),
        UnitSphericalPoint(0.30139703746681035, 0.8280815547305803, -0.47269521313549934)
    ]

    p1 = GI.Polygon([p1_pts])
    p2 = GI.Polygon([p2_pts])

    # This should not throw TracingError
    result = @test_nowarn GO.intersection(GO.Spherical(), p1, p2; target=GI.PolygonTrait())

    # Should produce non-empty intersection (polygons overlap)
    @test !isempty(result)

    # Verify the intersection has positive area
    total_area = sum(p -> GO.area(GO.Spherical(), p), result)
    @test total_area > 0
end