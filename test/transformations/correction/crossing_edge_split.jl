using Test
import GeoInterface as GI
import GeometryOps as GO
import GeometryOps: Spherical, RelateNG

# A ring that self-intersects on the sphere: an explicit figure-eight whose
# diagonals cross at (0, 0). `prepare` on the Spherical manifold rejects it;
# `CrossingEdgeSplit` is the documented remedy.
bowtie_pts = [(-10., -10.), (10., 10.), (10., -10.), (-10., 10.), (-10., -10.)]
bowtie = GI.Polygon([GI.LinearRing(bowtie_pts)])
clean_poly = GI.Polygon([GI.LinearRing([(0., 0.), (10., 0.), (10., 10.), (0., 10.), (0., 0.)])])

salg = RelateNG(; manifold = Spherical())

@testset "CrossingEdgeSplit" begin
    @testset "bowtie splits into two shells with even-odd containment" begin
        for pts in (bowtie_pts, reverse(bowtie_pts))
            poly = GI.Polygon([GI.LinearRing(pts)])
            @test_throws ArgumentError GO.prepare(salg, poly)
            fixed = GO.CrossingEdgeSplit()(poly)
            @test GI.trait(fixed) isa GI.MultiPolygonTrait
            @test GI.ngeom(fixed) == 2
            #-- each lobe is a triangle: the crossing vertex plus two originals
            @test sort([GI.npoint(g) for g in GI.getgeom(fixed)]) == [4, 4]
            #-- the constructed crossing vertex sits at the diagonals' crossing
            allpts = collect(GI.getpoint(GI.getexterior(first(GI.getgeom(fixed)))))
            @test any(p -> abs(GI.x(p)) < 1e-9 && abs(GI.y(p)) < 1e-9, allpts)
            #-- the repair validates clean
            @test GO.prepare(salg, fixed) isa GO.PreparedRelate
            #-- containment matches even-odd truth in both lobes and outside
            @test GO.intersects(salg, fixed, GI.Point(8., 1.))
            @test GO.intersects(salg, fixed, GI.Point(-8., -1.))
            @test !GO.intersects(salg, fixed, GI.Point(0., 5.))
            @test !GO.intersects(salg, fixed, GI.Point(0., -5.))
            @test !GO.intersects(salg, fixed, GI.Point(100., 40.))
        end
    end

    @testset "valid geometry is returned unchanged (identity, no copy)" begin
        @test GO.CrossingEdgeSplit()(clean_poly) === clean_poly
        mp = GI.MultiPolygon([clean_poly])
        @test GO.CrossingEdgeSplit()(mp) === mp
    end

    @testset "fix() integration" begin
        fixed = GO.fix(bowtie; corrections = [GO.CrossingEdgeSplit()])
        @test GI.trait(fixed) isa GI.MultiPolygonTrait
        @test GI.ngeom(fixed) == 2
    end

    @testset "holes are assigned to their enclosing shell loop" begin
        holed = GI.Polygon([
            GI.LinearRing(bowtie_pts),
            GI.LinearRing([(7., -1.), (9., -1.), (9., 1.), (7., 1.), (7., -1.)]),
        ])
        fixed = GO.CrossingEdgeSplit()(holed)
        @test GI.trait(fixed) isa GI.MultiPolygonTrait
        @test GI.ngeom(fixed) == 2
        @test sort([GI.nring(g) for g in GI.getgeom(fixed)]) == [1, 2]
        @test GO.prepare(salg, fixed) isa GO.PreparedRelate
        @test !GO.intersects(salg, fixed, GI.Point(8., 0.))     # inside the hole
        @test GO.intersects(salg, fixed, GI.Point(9.5, 0.))     # lobe, past the hole
        @test GO.intersects(salg, fixed, GI.Point(-8., -1.))    # the unholed lobe
    end

    @testset "a self-crossing hole splits into two hole loops" begin
        poly = GI.Polygon([
            GI.LinearRing([(-20., -20.), (20., -20.), (20., 20.), (-20., 20.), (-20., -20.)]),
            GI.LinearRing(bowtie_pts),
        ])
        fixed = GO.CrossingEdgeSplit()(poly)
        @test GI.trait(fixed) isa GI.PolygonTrait   # the shell did not split
        @test GI.nring(fixed) == 3
        @test GO.prepare(salg, fixed) isa GO.PreparedRelate
        @test !GO.intersects(salg, fixed, GI.Point(8., 1.))     # in a hole lobe
        @test !GO.intersects(salg, fixed, GI.Point(-8., -1.))   # in the other lobe
        @test GO.contains(salg, fixed, GI.Point(0., 5.))        # between the lobes
        @test GO.contains(salg, fixed, GI.Point(-15., 0.))      # ordinary interior
    end

    @testset "tangled multi-crossing topology throws instead of guessing" begin
        #-- a pentagram: every edge properly crosses two others
        star_angles = (90., 234., 18., 162., 306.)
        star = [(10cosd(a), 10sind(a)) for a in star_angles]
        pentagram = GI.Polygon([GI.LinearRing([star; [star[1]]])])
        @test_throws ArgumentError GO.CrossingEdgeSplit()(pentagram)
    end
end
