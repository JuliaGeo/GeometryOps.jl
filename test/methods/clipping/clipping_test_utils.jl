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
p27 = GI.Polygon([[[0.0, 0.0], [8.0, 0.0], [10.0, -1.0], [8.0, 1.0], [8.0, 2.0], [7.0, 5.0], [6.0, 4.0], [3.0, 5.0], [3.0, 3.0], [0.0, 0.0]]])
p28 = GI.Polygon([[[1.0, 1.0], [3.0, -1.0], [6.0, 2.0], [8.0, 0.0], [8.0, 4.0], [4.0, 4.0], [1.0, 1.0]]])
p29 = GI.Polygon([[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [3.0, 1.0], [1.0, 1.0], [0.0, 2.0], [0.0, 0.0]]])
p30 = GI.Polygon([[[4.0, 0.0], [3.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 2.0], [4.0, 2.0], [4.0, 0.0]]])
p31 = GI.Polygon([[[0.0, 0.0], [2.0, 1.0], [4.0, 0.0], [2.0, 4.0], [1.0, 2.0], [0.0, 3.0], [0.0, 0.0]]])
p32 = GI.Polygon([[[4.0, 3.0], [3.0, 2.0], [4.0, 2.0], [4.0, 3.0]]])
p33 = GI.Polygon([[[0.0, 0.0], [3.0, 0.0], [3.0, 3.0], [0.0, 3.0], [0.0, 0.0]]])
p34 = GI.Polygon([[[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0], [1.0, 0.0]]])
p35 = GI.Polygon([[[1.0, 0.0], [2.0, 0.0], [2.0, -1.0], [1.0, -1.0], [1.0, 0.0]]])
p36 = GI.Polygon([[[2.0, 1.0], [3.0, 0.0], [4.0, 1.0], [3.0, 3.0], [2.0, 1.0]]])
p37 = GI.Polygon([[[1.0, -1.0], [2.0, -1.0], [2.0, -2.0], [1.0, -2.0], [1.0, -1.0]]])
p38 = GI.Polygon([[(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)], [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]])
p39 = GI.Polygon([[(5.0, 0.0), (8.0, 0.0), (8.0, 3.0), (5.0, 3.0), (5.0, 0.0)], [(6.0, 1.0), (7.0, 1.0), (7.0, 2.0), (7.0, 1.0), (6.0, 1.0)]])
p40 = GI.Polygon([[(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0), (0.0, 0.0)], [(5.0, 5.0), (7.0, 5.0), (7.0, 7.0), (5.0, 7.0), (5.0, 5.0)]])
p41 = GI.Polygon([[(3.0, -1.0), (10.0, -1.0), (10.0, 9.0), (3.0, 9.0), (3.0, -1.0)], [(4.0, 3.0), (5.0, 3.0), (5.0, 4.0), (4.0, 4.0), (4.0, 3.0)], [(6.0, 1.0), (7.0, 1.0), (7.0, 2.0), (6.0, 2.0), (6.0, 1.0)]])
p42 = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)], [(1.0, 1.0), (3.0, 1.0), (3.0, 1.5), (1.0, 1.5), (1.0, 1.0)], [(1.0, 2.5), (3.0, 2.5), (3.0, 3.0), (1.0, 3.0), (1.0, 2.5)]])
p43 = GI.Polygon([[(2.0, -1.0), (5.0, -1.0), (5.0, 5.0), (2.0, 5.0), (2.0, -1.0)], [(3.5, 4.0), (4.5, 4.0), (4.5, 4.5), (3.5, 4.5), (3.5, 4.0)], [(3.5, 3.0), (4.5, 3.0), (4.5, 3.5), (3.5, 3.5), (3.5, 3.0)], [(3.5, 2.0), (4.5, 2.0), (4.5, 2.5), (3.5, 2.5), (3.5, 2.0)]])
p44 = GI.Polygon([[(2.0, -1.0), (5.0, -1.0), (5.0, 5.0), (2.0, 5.0), (2.0, -1.0)], [(3.5, 3.0), (4.5, 3.0), (4.5, 3.5), (3.5, 3.5), (3.5, 3.0)], [(3.5, 2.0), (4.5, 2.0), (4.5, 2.5), (3.5, 2.5), (3.5, 2.0)], [(3.5, 1.0), (4.5, 1.0), (4.5, 1.5), (3.5, 1.5), (3.5, 1.0)]])
p45 = GI.Polygon([[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0), (0.0, 0.0)]])
p46 = GI.Polygon([[(1.0, 1.0), (4.0, 1.0), (4.0, 4.0), (1.0, 4.0), (1.0, 1.0)], [(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0), (2.0, 2.0)]])



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
    (p17, p18, "p17", "p18", "First polygon with a hole (hole completly in second polygon), second without a hole"),
    (p18, p17, "p18", "p17", "First polygon with no hole, second with a hole (hole completly in first polygon)"),
    (p19, p20, "p19", "p20", "First polygon with a hole (hole not completly in second polygon), second without a hole"),
    (p42, p20, "p42", "p20", "First polygon with two holes (holes not completly in second polygon), second without a hole"),
    (p20, p19, "p20", "p19", "First polygon with no hole, second with a hole (hole not completly in first polygon)"),
    (p20, p42, "p20", "p42", "First polygon with no holes, second with two holes (holes not completly in second polygon)"),
    (p21, p22, "p21", "p22", "Polygons form cross, splitting each other"),
    (p23, p24, "p23", "p24", "Polygons are both donuts with intersecting holes"),
    (p25, p26, "p25", "p26", "Polygons both have two holes that intersect in various ways"),
    (p27, p28, "p27", "p28", "Figure 12 from Foster extension for degeneracies"),
    # (p29, p30, "p29", "p30", "Figure 13 from Foster extension for degeneracies"),
    (p31, p32, "p31", "p32", "Polygons touch at just one point"),
    (p33, p34, "p33", "p34", "One polygon inside of the other, sharing an edge"),
    (p33, p35, "p33", "p35", "Polygons outside of one another, sharing an edge"),
    (p33, p36, "p33", "p36", "All intersection points are V-intersections as defined by Foster"),
    (p33, p37, "p33", "p37", "Polygons are completly disjoint (no holes)"),
    (p38, p39, "p38", "p39", "Polygons are completly disjoint (both have one hole)"),
    (p40, p41, "p40", "p41", "Two overlapping polygons with three total holes in overlap region"),
    (p42, p43, "p42", "p43", "First polygon 2 holes, second polygon 3 holes. Holes do not overlap"),
    (p43, p42, "p43", "p42", "First polygon 3 holes, second polygon 2 holes. Holes do not overlap"),
    (p42, p44, "p42", "p43", "First polygon 2 holes, second polygon 3 holes. Holes do not overlap"),
    (p44, p42, "p43", "p42", "First polygon 3 holes, second polygon 2 holes. Holes do not overlap")

]
const ϵ = 1e-10
function compare_GO_LG_clipping(GO_f, LG_f, p1, p2)
    GO_result_list = GO_f(p1, p2; target = GI.PolygonTrait)
    LG_result_geom = LG_f(p1, p2)
    if LG_result_geom isa LG.GeometryCollection
        poly_list = LG.Polygon[]
        for g in GI.getgeom(LG_result_geom)
            g isa LG.Polygon && push!(poly_list, g)
        end
        LG_result_geom = LG.MultiPolygon(poly_list)
    end
    if isempty(GO_result_list) && (LG.isEmpty(LG_result_geom) || LG.area(LG_result_geom) == 0)
        return true
    end
    local GO_result_geom
    if length(GO_result_list)==1
        GO_result_geom = GO_result_list[1]
    else
        GO_result_geom = GI.MultiPolygon(GO_result_list)
    end
    diff_area = LG.area(LG.difference(GO_result_geom, LG_result_geom))
    return diff_area ≤ ϵ
end

function test_clipping(GO_f, LG_f, f_name)
    for (p1, p2, sg1, sg2, sdesc) in test_pairs
        println("$sg1 and $sg2")
        println(f_name)
        pass_test = compare_GO_LG_clipping(GO_f, LG_f, p1, p2)
        @test pass_test
        !pass_test && println("\n↑ TEST INFO: $sg1 $f_name $sg2 - $sdesc \n\n")
    end
end

@testset "Intersection" begin test_clipping(GO.intersection, LG.intersection, "intersection") end
@testset "Union" begin test_clipping(GO.union, LG.union, "union") end
@testset "Difference" begin test_clipping(GO.difference, LG.difference, "difference") end