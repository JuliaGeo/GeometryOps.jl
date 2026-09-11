using Test
import GeoInterface as GI
import GeometryOps as GO

# Wrap longitudes deliberately: planar bounding boxes and predicates give different
# answers here, so a missing manifold cannot pass these tests by accident.
function passthrough_box(x1, y1, x2, y2, ::Type{T} = Float64) where T
    pts = [(x1,y1), (x2,y1), (x2,y2), (x1,y2), (x1,y1)]
    GI.Polygon([[T.((mod(x + 180, 360) - 180, y)) for (x,y) in pts]])
end
passthrough_xyz(p) = GO.apply(GI.PointTrait(), p) do pt
    GO.UnitSpherical.UnitSphericalPoint(pt)
end
passthrough_area(m, polys) = sum(p -> GO.area(m, p), polys; init = 0.0)

@testset "Manifold and numeric type through clipping" begin
    for T in (Float32, Float64), xyz in (false, true)
        m = GO.Spherical()
        convertpoly = xyz ? passthrough_xyz : identity
        a = convertpoly(passthrough_box(170,-10,190,10,T))
        b = convertpoly(passthrough_box(175,-5,195,15,T))
        expected_type = xyz ? GO.UnitSpherical.UnitSphericalPoint{T} : Tuple{T,T}
        for op in (GO.intersection, GO.union, GO.difference)
            reference = op(m, a, b, T; target = GI.PolygonTrait())
            @test !isempty(reference)
            @test all(p -> all(pt -> pt isa expected_type, GI.getpoint(p)), reference)
            for aa in (a, GI.MultiPolygon([a,a])), bb in (b, GI.MultiPolygon([b,b]))
                result = op(m, aa, bb, T; target = GI.PolygonTrait())
                @test passthrough_area(m, result) ≈ passthrough_area(m, reference) rtol=(T === Float32 ? 2e-5 : 1e-12)
                @test all(p -> all(pt -> pt isa expected_type, GI.getpoint(p)), result)
                multi = op(m, aa, bb, T; target = GI.MultiPolygonTrait())
                @test GO.area(m, multi) ≈ passthrough_area(m, reference) rtol=(T === Float32 ? 2e-5 : 1e-12)
            end
        end
        for correction in (GO.UnionIntersectingPolygons(m,T), GO.DiffIntersectingPolygons(m,T))
            result = correction(GI.MultiPolygon([a,a]))
            @test GI.ngeom(result) == 1
            @test GO.area(m,result) ≈ GO.area(m,a)
            @test first(GI.getpoint(first(GI.getgeom(result)))) isa expected_type
        end
    end
end

@testset "Spherical holes" begin
    m = GO.Spherical()
    for xyz in (false,true)
        conv = xyz ? passthrough_xyz : identity
        outer = conv(passthrough_box(160,-20,200,20))
        hole = conv(passthrough_box(172,-8,188,8))
        island = conv(passthrough_box(176,-4,184,4))
        donut = GI.Polygon([GI.getexterior(outer),GI.getexterior(hole)])
        # Containment uses holes, including when every longitude box spans the dateline.
        @test isempty(GO.intersection(m,donut,island;target=GI.PolygonTrait()))
        united = GO.union(m,donut,island;target=GI.PolygonTrait())
        @test length(united) == 2
        @test passthrough_area(m,united) ≈ GO.area(m,donut)+GO.area(m,island)
        @test passthrough_area(m,GO.difference(m,outer,hole;target=GI.PolygonTrait())) ≈ GO.area(m,donut)
        # Overlapping holes exercise their union and intersection, and a partial hole
        # intersection creates a notch in the exterior rather than a closed hole.
        hole2 = conv(passthrough_box(180,-6,192,6))
        donut2 = GI.Polygon([GI.getexterior(outer),GI.getexterior(hole2)])
        for b in (donut2,conv(passthrough_box(182,-15,205,15)))
            i = GO.intersection(m,donut,b;target=GI.PolygonTrait())
            u = GO.union(m,donut,b;target=GI.PolygonTrait())
            d = GO.difference(m,donut,b;target=GI.PolygonTrait())
            @test passthrough_area(m,i)+passthrough_area(m,u) ≈ GO.area(m,donut)+GO.area(m,b) rtol=1e-10
            @test passthrough_area(m,i)+passthrough_area(m,d) ≈ GO.area(m,donut) rtol=1e-10
        end
    end
end

@testset "Spherical intersection points and accelerators" begin
    m = GO.Spherical()
    # The great-circle arc bulges north of both endpoints: vertex envelopes are disjoint.
    a = GI.Line([(-45.0,60.0),(45.0,60.0)])
    b = GI.Line([(0.0,65.0),(0.0,80.0)])
    for T in (Float32,Float64), xyz in (false,true)
        aa,bb = xyz ? (passthrough_xyz(a),passthrough_xyz(b)) : (a,b)
        result = GO.intersection_points(m,aa,bb,T)
        @test length(result) == 1
        @test first(result) isa (xyz ? GO.UnitSpherical.UnitSphericalPoint{T} : Tuple{T,T})
        @test result == GO.intersection_points(m,GO.AutoAccelerator(),aa,bb,T)
        @test result == GO.intersection(m,aa,bb,T;target=GI.PointTrait())
    end
    for accelerator in (GO.SingleSTRtree(),GO.SingleNaturalTree(),GO.DoubleNaturalTree(),GO.ThinnedDoubleNaturalTree())
        @test_throws ArgumentError GO.intersection_points(m,accelerator,a,b)
    end
end

@testset "Spherical cut preserves representation" begin
    for T in (Float32,Float64), xyz in (false,true)
        conv = xyz ? passthrough_xyz : identity
        p = conv(passthrough_box(170,-10,190,10,T))
        line = conv(GI.Line([T.((180,-20)),T.((180,20))]))
        pieces = GO.cut(GO.Spherical(),p,line,T)
        @test length(pieces) == 2
        @test passthrough_area(GO.Spherical(),pieces) ≈ GO.area(GO.Spherical(),p) rtol=(T === Float32 ? 2e-6 : 1e-12)
        @test all(poly -> all(pt -> pt isa (xyz ? GO.UnitSpherical.UnitSphericalPoint{T} : Tuple{T,T}), GI.getpoint(poly)),pieces)
    end
end

@testset "Correction metadata and empty multipolygons" begin
    crs = :test_crs
    p = passthrough_box(0,0,5,5)
    multi = GI.MultiPolygon([p,p]; crs)
    for correction in (GO.UnionIntersectingPolygons(),GO.DiffIntersectingPolygons())
        @test GI.crs(correction(multi)) == crs
    end
    empty_multi = GO.difference(p,p;target=GI.MultiPolygonTrait())
    @test GI.ngeom(empty_multi) == 0
    for op in (GO.union,GO.intersection,GO.difference)
        result = op(empty_multi,empty_multi;target=GI.MultiPolygonTrait())
        @test GI.ngeom(result) == 0
    end
end

include(joinpath(@__DIR__, "..", "..", "external", "s2geography", "s2geography.jl"))
using .S2Geog
using GeometryOpsTestHelpers: write_wkb

@testset "Spherical holes agree with S2" begin
    if s2_available()
        outer = passthrough_box(160,-20,200,20)
        hole = passthrough_box(172,-8,188,8)
        hole2 = passthrough_box(180,-6,192,6)
        donut = GI.Polygon([GI.getexterior(outer),GI.getexterior(hole)])
        donut2 = GI.Polygon([GI.getexterior(outer),GI.getexterior(hole2)])
        m = GO.Spherical(radius = S2Geog.S2_RADIUS)
        for b in (donut2,passthrough_box(182,-15,205,15),passthrough_box(176,-4,184,4)),
            (name,op) in ((:intersection,GO.intersection),(:union,GO.union),(:difference,GO.difference))
            expected = s2_area(s2_overlay(name,write_wkb(donut),write_wkb(b)))
            scale = s2_area(write_wkb(donut)) + s2_area(write_wkb(b))
            for conv in (identity,passthrough_xyz)
                result = op(m,conv(donut),conv(b);target=GI.PolygonTrait())
                @test abs(passthrough_area(m,result)-expected) ≤ 1e-10*scale
            end
        end
    else
        @test_skip false
    end
end

@testset "Cut line-order pairing on concave polygons" begin
    points = [(0.,0.),(6.,0.),(6.,6.),(4.,6.),(4.,2.),(2.,2.),(2.,6.),(0.,6.)]
    for m in (GO.Planar(), GO.Spherical()), reverse_ring in (false,true), offset in eachindex(points)
        pts = circshift(reverse_ring ? reverse(points) : points, offset)
        p = GI.Polygon([vcat(pts,[pts[1]])])
        fully_crossing = GI.Line([(-1.,4.),(7.,4.)])
        pieces = GO.cut(m,p,fully_crossing)
        @test length(pieces) == 3
        @test passthrough_area(m,pieces) ≈ GO.area(m,p) rtol=1e-12
        # Both endpoints lie in different lobes. Two crossings do not mean a full cut.
        partial = GI.Line([(1.,4.),(5.,4.)])
        @test GI.coordinates(only(GO.cut(m,p,partial))) == GI.coordinates(p)
    end
end

@testset "Mixed spherical representations" begin
    a = passthrough_box(170,-10,190,10)
    b = passthrough_box(175,-5,195,15)
    m = GO.Spherical()
    for xyz_a in (false,true), xyz_b in (false,true), op in (GO.intersection,GO.union,GO.difference)
        aa = xyz_a ? passthrough_xyz(a) : a
        bb = xyz_b ? passthrough_xyz(b) : b
        expected_type = xyz_a ? GO.UnitSpherical.UnitSphericalPoint{Float64} : Tuple{Float64,Float64}
        reference = op(m,a,b;target=GI.PolygonTrait())
        for firstarg in (aa,GI.MultiPolygon([aa])), lastarg in (bb,GI.MultiPolygon([bb]))
            result = op(m,firstarg,lastarg;target=GI.PolygonTrait())
            @test passthrough_area(m,result) ≈ passthrough_area(m,reference) rtol=1e-12
            @test all(p -> all(pt -> pt isa expected_type,GI.getpoint(p)),result)
            multi = op(m,firstarg,lastarg;target=GI.MultiPolygonTrait())
            @test GO.area(m,multi) ≈ passthrough_area(m,reference) rtol=1e-12
            @test all(p -> all(pt -> pt isa expected_type,GI.getpoint(p)),GI.getpolygon(multi))
        end
    end
end

@testset "Planar numeric-type forwarding" begin
    a = passthrough_box(0,0,2,2)
    b = passthrough_box(1,-1,3,1)
    for op in (GO.union,GO.intersection,GO.difference), prefix in ((),(GO.Planar(),))
        ps = op(prefix...,a,b,Float32;target=GI.PolygonTrait())
        @test first(GI.getpoint(first(ps))) isa Tuple{Float32,Float32}
    end
end

@testset "Mixed representation passthrough is bit-exact" begin
    m = GO.Spherical()
    a = GI.Polygon([[ (170.123456789,-2.123456789), (174.987654321,-1.987654321),
        (173.765432109,3.123456789), (170.123456789,-2.123456789) ]])
    outer = passthrough_xyz(passthrough_box(165,-10,180,10))
    disjoint = passthrough_xyz(passthrough_box(140,-5,145,5))
    ma = GI.MultiPolygon([a])
    result = GO.intersection(m,ma,outer;target=GI.PolygonTrait())
    @test collect(GI.getpoint(only(result))) == collect(GI.getpoint(a))
    for result in (GO.union(m,ma,disjoint;target=GI.PolygonTrait()), GO.difference(m,ma,disjoint;target=GI.PolygonTrait()))
        pts = collect(Iterators.flatten(GI.getpoint.(result)))
        @test all(p -> p in pts, GI.getpoint(a))
    end
end
