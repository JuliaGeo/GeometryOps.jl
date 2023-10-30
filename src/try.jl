import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

using GeoInterface.Extents: Extents
using GeoInterface

using LinearAlgebra
import ExactPredicates
using Random

const TuplePoint = Tuple{Float64,Float64}
const Edge = Tuple{TuplePoint,TuplePoint}

function compareGO_LG(poly_1,poly_2, ϵ)
    p1 = GI.Polygon([poly_1])
    p2 = GI.Polygon([poly_2])
    LG_p1p2 = LG.intersection(p1, p2)
    GO_p1p2 = GO.intersection(p1, p2)
    # if I do GI.equals(LG_p1p2, GO_p1p2) here I get false even when it is true
    # if I wrap output of LG.intersection in GI Polygon, it still gets false
    # the only thing that returns true is turning both in LG polygon
    if length(GO_p1p2)==1
        inter_GO = LG.Polygon(convert_tuple_to_array(GO_p1p2))
    else
        temp = convert_tuple_to_array(GO_p1p2)
        inter_GO = LG.MultiPolygon([temp])
    end
    return LG.area(LG.difference(inter_GO, LG_p1p2)) < ϵ
end

function convert_tuple_to_array(tups)
    return_polys = Array{Array{Array{Float64, 1}, 1}, 1}(undef, 0)
    for polygon in tups
        pt_list = Array{Array{Float64, 1}, 1}(undef, 0)
        for point in polygon
            push!(pt_list, [point[1], point[2]])
        end
        push!(return_polys, pt_list)
    end
    return return_polys
end

p1 = [[0.0, 0.0], [5.0, 5.0], [10.0, 0.0], [5.0, -5.0], [0.0, 0.0]]
p2 = [[3.0, 0.0], [8.0, 5.0], [13.0, 0.0], [8.0, -5.0], [3.0, 0.0]]

poly_1 = [(4.526700198111509, 3.4853728532584696), (2.630732683726619, -4.126134282323841),
     (-0.7638522032421201, -4.418734350277446), (-4.367920073785058, -0.2962672719707883),
     (4.526700198111509, 3.4853728532584696)]

poly_2 = [(5.895141140952208, -0.8095078714426418), (2.8634927670695283, -4.625511746720306),
     (-1.790623183259246, -4.138092164660989), (-3.9856656502985843, -0.5275687876429914),
     (-2.554809853598822, 3.553455552936806), (1.1865909598835922, 4.984203644564732),
     (5.895141140952208, -0.8095078714426418)]

# example polygon from greiner paper
p3 = [(0.0, 0.0), (0.0, 4.0), (7.0, 4.0), (7.0, 0.0), (0.0, 0.0)]
p4 = [(1.0, -3.0), (1.0, 1.0), (3.5, -1.5), (6.0, 1.0), (6.0, -3.0), (1.0, -3.0)]

# polygons made with high spikiness so concave
p5 = [(1.2938349167338743, -3.175128530227131), (-2.073885870841754, -1.6247711001754137),
(-5.787437985975053, 0.06570713422599561), (-2.1308128111898093, 5.426689675486368),
(2.3058074184797244, 6.926652158268195), (1.2938349167338743, -3.175128530227131)]
p6 = [(-2.1902469793743924, -1.9576242117579579), (-4.726006206053999, 1.3907098941556428),
(-3.165301985923147, 2.847612825874245), (-2.5529280962099428, 4.395492123980911),
(0.5677700216973937, 6.344638314896882), (3.982554842356183, 4.853519613487035),
(5.251193948893394, 0.9343031382106848), (5.53045582244555, -3.0101433691361734),
(-2.1902469793743924, -1.9576242117579579)]


# compareGO_LG(p5, p6, 1e-5)
poly_a = GI.Polygon([p5])
poly_b = GI.Polygon([p6])
include(raw"C:\Users\Lana\.julia\dev\GeometryOps\src\utils.jl")
edges_b = to_edges(poly_b);

include(raw"C:\Users\Lana\.julia\dev\GeometryOps\src\methods\bools.jl")

point_in_polygon(edges_b[1][1], poly_a)
# GO.intersection(poly_a, poly_b)

# p2 = LG.Point([0.0, 1.0])
# mp3 = LG.MultiPoint([p2])
# GO.equals(p2, mp3)
