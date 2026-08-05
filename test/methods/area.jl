using Test
import GeoInterface as GI
import GeometryOps as GO 
import LibGEOS as LG
import GeoFormatTypes
import Proj
using GeometryOpsTestHelpers

pt = LG.Point([0.0, 0.0])
empty_pt = LG.readgeom("POINT EMPTY")
mpt = LG.MultiPoint([[0.0, 0.0], [1.0, 0.0]])
empty_mpt = LG.readgeom("MULTIPOINT EMPTY")
l1 = LG.LineString([[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]])
empty_l = LG.readgeom("LINESTRING EMPTY")
ml1 = LG.MultiLineString([[[0.0, 0.0], [0.5, 0.5], [1.0, 0.5]], [[0.0, 0.0], [0.0, 0.1]]])
empty_ml = LG.readgeom("MULTILINESTRING EMPTY")
empty_l = LG.readgeom("LINESTRING EMPTY")
r1 = LG.LinearRing([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 0.0]])
empty_r = LG.readgeom("LINEARRING EMPTY")
p1 = LG.Polygon([
    [[10.0, 0.0], [30.0, 0.0], [30.0, 20.0], [10.0, 20.0], [10.0, 0.0]],
])
p2 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [15.0, 11.0], [25.0, 11.0], [25.0, 1.0], [15.0, 1.0]],
])
p3 = LG.Polygon([
    [[10.0, 0.0], [10.0, 20.0], [30.0, 20.0], [30.0, 0.0], [10.0, 0.0]],
    [[15.0, 1.0], [25.0, 1.0], [25.0, 11.0], [15.0, 11.0], [15.0, 1.0]],
])
p4 = LG.Polygon([
    [
        [0.0, 5.0], [2.0, 2.0], [5.0, 2.0], [2.0, -2.0], [5.0, -5.0],
        [0.0, -2.0], [-5.0, -5.0], [-2.0, -2.0], [-5.0, 2.0], [-2.0, 2.0],
        [0.0, 5.0],
    ],
])
empty_p = LG.readgeom("POLYGON EMPTY")
mp1 = LG.MultiPolygon([p2, p4])
empty_mp = LG.readgeom("MULTIPOLYGON EMPTY")
c = LG.GeometryCollection([p1, p2, r1])
c_with_epty_l = LG.GeometryCollection([p1, p2, r1, empty_l])
empty_c = LG.readgeom("GEOMETRYCOLLECTION EMPTY")

@testset_implementations "That handle empty geoms" begin 
    @test GO.area($empty_pt) == LG.area($empty_pt) == 0
    @test GO.area($empty_mpt) == LG.area($empty_mpt) == 0
    @test GO.area($empty_l) == LG.area($empty_l) == 0
    @test GO.area($empty_ml) == LG.area($empty_ml) == 0
    @test GO.area($empty_r) == LG.area($empty_r) == 0
    # Empty polygon
    @test GO.signed_area($empty_p) == 0
    @test GO.area($empty_p) == LG.area($empty_p) == 0
    # Empty multipolygon
    @test GO.area($empty_mp) == LG.area($empty_mp) == 0
    # Empty collection
    @test GO.area(c_with_epty_l) == LG.area(c_with_epty_l)
    @test GO.area(c_with_epty_l, Float32) isa Float32
    @test GO.area(empty_c) == LG.area(empty_c) == 0
end

@testset "With GeometryCollection" begin 
    # Geometry collection summed area
    @test GO.area(c) == LG.area(c)
    @test GO.area(c, Float32) isa Float32
end  

@testset_implementations "all" begin 
    # Points, lines, and rings have zero area
    @test GO.area($pt) == GO.signed_area($pt) == LG.area($pt) == 0
    @test GO.area($pt) isa Float64
    @test GO.signed_area($pt, Float32) isa Float32
    @test GO.signed_area($pt) isa Float64
    @test GO.area($pt, Float32) isa Float32
    @test GO.area($mpt) == GO.signed_area($mpt) == LG.area($mpt) == 0
    @test GO.area($l1) == GO.signed_area($l1) == LG.area($l1) == 0
    @test GO.area($ml1) == GO.signed_area($ml1) == LG.area($ml1) == 0
    @test GO.area($r1) == GO.signed_area($r1) == LG.area($r1) == 0

    # Polygons have non-zero area
    # CCW polygons have positive signed area
    @test GO.area($p1) == GO.signed_area($p1) == LG.area($p1)
    @test GO.signed_area($p1) > 0
    # Float32 calculations
    @test GO.area($p1) isa Float64
    @test GO.area($p1, Float32) isa Float32
    # CW polygons have negative signed area
    a2 = LG.area($p2)
    @test GO.area($p2) == a2
    @test GO.signed_area($p2) == -a2
    # Winding order of holes doesn't affect sign of signed area
    @test GO.signed_area($p3) == -a2
    # Concave polygon correctly calculates area
    a4 = LG.area($p4)
    @test GO.area($p4) == a4
    @test GO.signed_area($p4) == -a4

    # Multipolygon calculations work
    @test GO.area($mp1) == a2 + a4
    @test GO.area($mp1, Float32) isa Float32
end


highlat_poly = LG.Polygon([[[70., 70.], [70., 80.], [80., 80.], [80., 70.], [70., 70.]]])

@testset "CRS-aware automatic area" begin
    geographic_square = GI.Polygon(
        [[
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        ]];
        crs=GeoFormatTypes.EPSG(4326),
    )
    @test GI.crstrait(geographic_square) isa GI.UnknownTrait
    @test GO.area(geographic_square) == GO.area(GO.Geodesic(), geographic_square)
    @test GO.area(geographic_square) ≈ 1.2308778361469452e10 rtol=1e-10
    @test GO.area(geographic_square) != GO.area(GO.Spherical(), geographic_square)

    foot_side = 3.280833333333333
    projected_square = GI.Polygon(
        [[
            (0.0, 0.0), (foot_side, 0.0), (foot_side, foot_side),
            (0.0, foot_side), (0.0, 0.0),
        ]];
        crs=GeoFormatTypes.EPSG(2263),
    )
    @test GO.area(projected_square) == GO.area(GO.Planar(), projected_square) == foot_side^2

    threaded_crsless_collection = GI.GeometryCollection(fill(projected_square, 256))
    @test GO.area(threaded_crsless_collection; threaded=true) ≈ 256 * foot_side^2

    @test GO.area([geographic_square]) == GO.area(GO.Planar(), geographic_square)
    @test GO.area([geographic_square, projected_square]) ==
          GO.area(GO.Planar(), [geographic_square, projected_square])
    crsless_feature_collection = GI.FeatureCollection([
        GI.Feature(geographic_square),
        GI.Feature(projected_square),
    ])
    @test GO.area(crsless_feature_collection) ==
          GO.area(GO.Planar(), crsless_feature_collection)

    plain_square = GI.Polygon([[
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
    ]])
    geographic_feature = GI.Feature(plain_square; crs=GeoFormatTypes.EPSG(4326))
    crsless_feature_collection = GI.FeatureCollection([geographic_feature])
    @test GO.area(crsless_feature_collection) ==
          GO.area(GO.Planar(), crsless_feature_collection)

    geographic_feature_collection = GI.FeatureCollection(
        [GI.Feature(plain_square)]; crs=GeoFormatTypes.EPSG(4326),
    )
    @test GO.area(geographic_feature_collection) ==
          GO.area(GO.Geodesic(), geographic_feature_collection)

    grads_square = GI.Polygon([[
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
    ]]; crs=GeoFormatTypes.EPSG(4807))
    @test GO.area(grads_square) ≈ 9.969411479067013e9 rtol=1e-10

    proj_crs_projected_square = GI.Polygon(
        [[
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        ]];
        crs=Proj.CRS(GeoFormatTypes.EPSG(2263)),
    )
    @test GO.area(proj_crs_projected_square) ==
          GO.area(GO.Planar(), proj_crs_projected_square) == 1.0

    mixed_crs_collection = GI.GeometryCollection([geographic_square, proj_crs_projected_square])
    @test isnothing(GI.crs(mixed_crs_collection))
    @test GO.area(mixed_crs_collection) == GO.area(GO.Planar(), mixed_crs_collection)

    crsless_multipolygon = GI.MultiPolygon([geographic_square, plain_square])
    @test isnothing(GI.crs(crsless_multipolygon))
    @test GO.area(crsless_multipolygon) == GO.area(GO.Planar(), crsless_multipolygon)

    geographic_collection = GI.GeometryCollection(
        [plain_square, projected_square]; crs=GeoFormatTypes.EPSG(4326),
    )
    @test GO.area(geographic_collection) ==
          GO.area(GO.Geodesic(), geographic_collection)

    spherical_crs_square = GI.Polygon(
        [[
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0),
        ]];
        crs=GeoFormatTypes.ProjString("+proj=longlat +R=6371000 +type=crs"),
    )
    spherical_crs_area = GO.area(spherical_crs_square)
    @test isfinite(spherical_crs_area)
    @test spherical_crs_area ≈ 1.2363997753679907e10 rtol=1e-10

    geographic_point = GI.Point((0.0, 0.0); crs=GeoFormatTypes.EPSG(4326))
    geographic_multipoint = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)]; crs=GeoFormatTypes.EPSG(4326))
    geographic_linestring = GI.LineString([(0.0, 0.0), (1.0, 1.0)]; crs=GeoFormatTypes.EPSG(4326))
    geographic_ring = GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)]; crs=GeoFormatTypes.EPSG(4326))
    @test GO.area(geographic_point) == 0.0
    @test GO.area(geographic_multipoint) == 0.0
    @test GO.area(geographic_linestring) == 0.0
    @test GO.area(geographic_ring) == 0.0
end

@testset "Geodesic polygon area" begin
    crs = GeoFormatTypes.EPSG(4326)
    exterior = [
        (0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0),
    ]
    hole = [
        (0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5),
    ]
    exterior_poly = GI.Polygon([exterior]; crs)
    reversed_exterior_poly = GI.Polygon([reverse(exterior)]; crs)
    hole_poly = GI.Polygon([hole]; crs)
    poly_with_hole = GI.Polygon([exterior, hole]; crs)
    poly_with_reversed_hole = GI.Polygon([exterior, reverse(hole)]; crs)

    exterior_area = GO.area(GO.Geodesic(), exterior_poly)
    hole_area = GO.area(GO.Geodesic(), hole_poly)
    @test GO.area(GO.Geodesic(), exterior_poly; init=2.0) == exterior_area + 2.0
    @test GO.area(GO.Geodesic(), reversed_exterior_poly) ≈ exterior_area rtol=1e-10
    @test GO.area(GO.Geodesic(), poly_with_hole) ≈ exterior_area - hole_area rtol=1e-10
    @test GO.area(GO.Geodesic(), poly_with_reversed_hole) ≈ exterior_area - hole_area rtol=1e-10

    second_poly = GI.Polygon([[
        (3.0, 0.0), (4.0, 0.0), (4.0, 1.0), (3.0, 1.0), (3.0, 0.0),
    ]]; crs)
    multipoly = GI.MultiPolygon([exterior_poly, second_poly]; crs)
    @test GO.area(GO.Geodesic(), multipoly) ≈
          exterior_area + GO.area(GO.Geodesic(), second_poly) rtol=1e-10

    threaded_multipoly = GI.MultiPolygon(fill(second_poly, 256); crs)
    threaded_area = GO.area(GO.Geodesic(), threaded_multipoly; threaded=true)
    @test GO.area(GO.Geodesic(), threaded_multipoly; threaded=true, init=2.0) ==
          threaded_area + 2.0
end

@testset_implementations "Spherical/geodesic" begin
    @test GO.area(GO.Planar(), $highlat_poly) == 100
    @test GO.area(GO.Planar(), $highlat_poly) < GO.area(GO.Geodesic(), $highlat_poly)
end

@testset "GirardSphericalArea algorithm type" begin
    # Test that the algorithm type exists and is a SingleManifoldAlgorithm for Spherical
    @test GO.NaiveTriangulatedSphericalArea <: GO.GeometryOpsCore.SingleManifoldAlgorithm{<:GO.Spherical}

    # Test construction with default manifold
    alg = GO.NaiveTriangulatedSphericalArea()
    @test alg isa GO.NaiveTriangulatedSphericalArea

    # Test manifold accessor
    @test GO.GeometryOpsCore.manifold(alg) isa GO.Spherical
end

@testset "Spherical triangle area" begin
    using GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic

    # Test a known spherical triangle: 1/8 of the unit sphere
    # Vertices at (0,0), (90,0), (0,90) in lon,lat
    # This is an octant of the sphere with area = 4πR²/8 = πR²/2
    p1 = UnitSphereFromGeographic()((0.0, 0.0))   # lon=0, lat=0
    p2 = UnitSphereFromGeographic()((90.0, 0.0))  # lon=90, lat=0
    p3 = UnitSphereFromGeographic()((0.0, 90.0))  # lon=0, lat=90 (north pole)

    # On unit sphere, area should be π/2
    area = GO._spherical_triangle_area(GO.Eriksson(), p1, p2, p3)
    @test area ≈ π/2 atol=1e-10

    # Test with reversed winding (should give same absolute area)
    area_rev = GO._spherical_triangle_area(GO.Eriksson(), p1, p3, p2)
    @test abs(area_rev) ≈ π/2 atol=1e-10
end

@testset "Spherical polygon area" begin
    unit_sphere = GO.Spherical(radius=1.0)

    # Simple polygon: octant (1/8 of sphere)
    # On unit sphere, area should be π/2
    octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])
    @test GO.area(unit_sphere, octant) ≈ π/2 atol=1e-10

    # Test a small polygon where spherical area should be positive
    small_poly = GI.Polygon([[
        (0.0, 0.0), (0.001, 0.0), (0.001, 0.001), (0.0, 0.001), (0.0, 0.0)
    ]])
    @test GO.area(unit_sphere, small_poly) > 0
end

@testset "Spherical polygon area with holes" begin
    unit_sphere = GO.Spherical(radius=1.0)

    # Polygon with a hole - the hole should subtract from the area
    exterior = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]
    hole = [(2.0, 2.0), (2.0, 3.0), (3.0, 3.0), (3.0, 2.0), (2.0, 2.0)]
    poly_with_hole = GI.Polygon([exterior, hole])
    poly_no_hole = GI.Polygon([exterior])
    hole_poly = GI.Polygon([hole])

    area_with_hole = GO.area(unit_sphere, poly_with_hole)
    area_exterior = GO.area(unit_sphere, poly_no_hole)
    area_hole = GO.area(unit_sphere, hole_poly)

    # Area with hole should be exterior minus hole
    @test area_with_hole ≈ area_exterior - area_hole atol=1e-10
end

@testset "area(Spherical(), geom) dispatch" begin
    # Test that Spherical manifold dispatches correctly
    octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])

    # Default Earth radius from Spherical()
    spherical = GO.Spherical()
    R = spherical.radius

    # Area should be (π/2) * R²
    expected_area = (π/2) * R^2
    computed_area = GO.area(spherical, octant)

    @test computed_area ≈ expected_area rtol=1e-10

    # Test with custom radius
    custom_spherical = GO.Spherical(radius=1.0)
    unit_area = GO.area(custom_spherical, octant)
    @test unit_area ≈ π/2 rtol=1e-10

    # Test that points/lines return zero
    pt = GI.Point((0.0, 0.0))
    @test GO.area(spherical, pt) == 0.0

    line = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
    @test GO.area(spherical, line) == 0.0
end

@testset "Spherical area integration tests" begin
    @testset "Known spherical areas" begin
        # Test 1: Octant (1/8 of sphere) = 4πR²/8 = πR²/2
        octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])
        R = GO.Spherical().radius
        @test GO.area(GO.Spherical(), octant) ≈ (π/2) * R^2 rtol=1e-8

        # Test 2: Very small polygon (should be positive and small)
        # At equator, 1° ≈ 111km, so 0.01° ≈ 1.1km
        tiny = GI.Polygon([[
            (0.0, 0.0), (0.01, 0.0), (0.01, 0.01), (0.0, 0.01), (0.0, 0.0)
        ]])
        tiny_spherical = GO.area(GO.Spherical(radius=1.0), tiny)
        # In radians: 0.01° = 0.01 * π/180 ≈ 1.745e-4 rad
        # Planar approx: (1.745e-4)² ≈ 3e-8 on unit sphere
        @test tiny_spherical > 0
        @test tiny_spherical < 1e-6  # Should be very small
    end

    @testset "Spherical vs Planar comparison" begin
        # High latitude polygon should have different spherical vs planar area
        highlat = GI.Polygon([[(70.0, 70.0), (80.0, 70.0), (80.0, 80.0), (70.0, 80.0), (70.0, 70.0)]])

        planar_area = GO.area(GO.Planar(), highlat)
        spherical_area = GO.area(GO.Spherical(), highlat)

        # Spherical area should be different from planar (in degrees²)
        # Spherical returns m², planar returns degrees²
        @test planar_area != spherical_area
    end

    @testset "MultiPolygon spherical area" begin
        poly1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
        poly2 = GI.Polygon([[(10.0, 10.0), (11.0, 10.0), (11.0, 11.0), (10.0, 11.0), (10.0, 10.0)]])
        mpoly = GI.MultiPolygon([poly1, poly2])

        area1 = GO.area(GO.Spherical(), poly1)
        area2 = GO.area(GO.Spherical(), poly2)
        multi_area = GO.area(GO.Spherical(), mpoly)

        @test multi_area ≈ area1 + area2 rtol=1e-10
    end

    @testset "Empty geometry spherical area" begin
        # Use LibGEOS to create empty polygon (GeoInterface.Polygon doesn't support empty rings)
        empty_poly = LG.readgeom("POLYGON EMPTY")
        @test GO.area(GO.Spherical(), empty_poly) == 0.0
    end

    @testset "Custom radius" begin
        # Mars radius approximately
        mars_radius = 3389.5e3  # meters
        mars = GO.Spherical(radius=mars_radius)

        octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])
        mars_area = GO.area(mars, octant)

        @test mars_area ≈ (π/2) * mars_radius^2 rtol=1e-8
    end

    @testset "area(NaiveTriangulatedSphericalArea(), geom) direct call" begin
        octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])

        # Direct algorithm call should use default Spherical manifold
        alg = GO.NaiveTriangulatedSphericalArea()
        R = GO.Spherical().radius

        expected_area = (π/2) * R^2
        computed_area = GO.area(alg, octant)

        @test computed_area ≈ expected_area rtol=1e-10
    end

    @testset "Type stability" begin
        # https://github.com/JuliaGeo/GeometryOps.jl/issues/407
        using GeometryOps.UnitSpherical: UnitSphericalPoint
        octant = GI.Polygon([[(0.0, 0.0), (90.0, 0.0), (0.0, 90.0), (0.0, 0.0)]])
        usp_octant = GI.Polygon([GI.LinearRing([
            UnitSphericalPoint(0.0, 0.0, 1.0), UnitSphericalPoint(1.0, 0.0, 0.0),
            UnitSphericalPoint(0.0, 1.0, 0.0), UnitSphericalPoint(0.0, 0.0, 1.0),
        ])])
        R = GO.Spherical().radius
        expected_area = (π/2) * R^2

        @test only(Base.return_types(GO.area, (GO.Spherical{Float64}, typeof(octant)))) == Float64
        @test only(Base.return_types(GO.area, (GO.Spherical{Float64}, typeof(usp_octant)))) == Float64
        @test only(Base.return_types(GO.area, (GO.NaiveTriangulatedSphericalArea{GO.Spherical{Float64}, GO.Eriksson}, typeof(octant)))) == Float64
        @test only(Base.return_types(GO.area, (GO.Spherical{Float64}, typeof(octant), Type{Float32}))) == Float32

        @test (@inferred GO.area(GO.Spherical(), octant)) ≈ expected_area rtol=1e-10
        @test (@inferred GO.area(GO.Spherical(), usp_octant)) ≈ expected_area rtol=1e-10
    end
end