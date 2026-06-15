using Test
import GeoInterface as GI
import GeometryOps as GO
import GeometryOps: Spherical, RelateNG

# The edge (0,0)→(180,0) maps to the antipodal unit vectors (1,0,0) and
# (-1,0,0), which the spherical kernel rejects.
antipodal_poly = GI.Polygon([GI.LinearRing([(0., 0.), (180., 0.), (90., 80.), (0., 0.)])])
clean_poly = GI.Polygon([GI.LinearRing([(0., 0.), (10., 0.), (10., 10.), (0., 10.), (0., 0.)])])

@testset "AntipodalEdgeSplit" begin
    alg = RelateNG(; manifold = Spherical())

    #-- without correction the spherical kernel refuses the antipodal edge
    @test_throws ArgumentError GO.relate(alg, antipodal_poly, GI.Point(10., 10.))

    #-- the correction inserts the lon/lat midpoint (90, 0): one extra vertex
    fixed = GO.AntipodalEdgeSplit()(antipodal_poly)
    @test GI.npoint(fixed) == GI.npoint(antipodal_poly) + 1
    pts = [(GI.x(p), GI.y(p)) for p in GI.getpoint(GI.getexterior(fixed))]
    @test (90.0, 0.0) in pts

    #-- after correction the relate runs and is correct (the big triangle
    #-- north of the equator contains the near-equator point)
    @test GO.relate(alg, fixed, GI.Point(10., 10.), "T*****FF*")

    #-- also reachable through `fix`
    @test GI.npoint(GO.fix(antipodal_poly; corrections = [GO.AntipodalEdgeSplit()])) ==
        GI.npoint(antipodal_poly) + 1

    #-- a geometry with no antipodal edge keeps its vertices
    @test GI.npoint(GO.AntipodalEdgeSplit()(clean_poly)) == GI.npoint(clean_poly)
end
