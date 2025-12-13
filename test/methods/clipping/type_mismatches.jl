using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import ArchGDAL as AG
using ..TestHelpers


p1 = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])])
p1_bothcrs = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]; crs = 1)], crs = 1)
p1_topcrs = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)])], crs = 1)
p1_extent = GO.tuples(p1_bothcrs; calc_extent = true)
p1_bothcrs_extent = GO.tuples(p1_bothcrs; calc_extent = true)
p1_multi = GI.MultiPolygon([p1])
p1_multi_bothcrs = GI.MultiPolygon([p1_bothcrs]; crs = 1)
p1_multi_topcrs = GI.MultiPolygon([p1_topcrs]; crs = 1)
p1_multi_extent = GO.tuples(p1_multi_bothcrs; calc_extent = true)
p1_multi_bothcrs_extent = GO.tuples(p1_multi_bothcrs; calc_extent = true)
p1_multi_topcrs_extent = GO.tuples(p1_multi_topcrs; calc_extent = true)

p2 = GO.transform(x -> x .+ (0.5, 0.0), p1)
p2_bothcrs = GO.transform(x -> x .+ (0.5, 0.0), p1_bothcrs)
p2_topcrs = GO.transform(x -> x .+ (0.5, 0.0), p1_topcrs)
p2_extent = GO.tuples(p2_bothcrs; calc_extent = true)
p2_bothcrs_extent = GO.tuples(p2_bothcrs; calc_extent = true)
p2_multi = GI.MultiPolygon([p2])
p2_multi_bothcrs = GI.MultiPolygon([p2_bothcrs]; crs = 1)
p2_multi_topcrs = GI.MultiPolygon([p2_topcrs]; crs = 1)
p2_multi_extent = GO.tuples(p2_multi_bothcrs; calc_extent = true)
p2_multi_bothcrs_extent = GO.tuples(p2_multi_bothcrs; calc_extent = true)
p2_multi_topcrs_extent = GO.tuples(p2_multi_topcrs; calc_extent = true)

p1s = zip(["p1", "p1_bothcrs", "p1_topcrs", "p1_extent", "p1_multi", "p1_multi_bothcrs", "p1_multi_topcrs", "p1_multi_extent", "p1_multi_bothcrs_extent", "p1_multi_topcrs_extent"], [p1, p1_bothcrs, p1_topcrs, p1_extent, p1_multi, p1_multi_bothcrs, p1_multi_topcrs, p1_multi_extent, p1_multi_bothcrs_extent, p1_multi_topcrs_extent])
p2s = zip(["p2", "p2_bothcrs", "p2_topcrs", "p2_extent", "p2_multi", "p2_multi_bothcrs", "p2_multi_topcrs", "p2_multi_extent", "p2_multi_bothcrs_extent", "p2_multi_topcrs_extent"], [p2, p2_bothcrs, p2_topcrs, p2_extent, p2_multi, p2_multi_bothcrs, p2_multi_topcrs, p2_multi_extent, p2_multi_bothcrs_extent, p2_multi_topcrs_extent])

function _test_falseiferror(f, args...; kwargs...)
    try
        f(args...; kwargs...)
        return true
    catch e
        return false
    end
end

@testset "Type mismatches" begin
    for (fname, func) in zip(["intersection", "union", "difference"], [GO.intersection, GO.union, GO.difference])
        @testset "$fname" begin
            for ((p1_name, p1_geom), (p2_name, p2_geom)) in Iterators.product(p1s, p2s)
                @testset_implementations "$fname $p1_name x $p2_name" begin
                    result = _test_falseiferror(func, $p1_geom, $p2_geom; target = GI.PolygonTrait)
                    @test result
                end
            end
        end
    end
end

@testset "Specifically ArchGDAL vs other things" begin

    # Create simple test polygons as GI.Polygon
    coords = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    gi_poly1 = GI.Polygon([GI.LinearRing(coords)])

    # Create a second overlapping polygon
    coords2 = [(0.5, 0.0), (1.5, 0.0), (1.5, 1.0), (0.5, 1.0), (0.5, 0.0)]
    gi_poly2 = GI.Polygon([GI.LinearRing(coords2)])

    # Create a polygon with CRS
    gi_poly_crs = GI.Polygon([GI.LinearRing(coords)]; crs=4326)

    # Create MultiPolygons
    gi_multi1 = GI.MultiPolygon([gi_poly1])
    gi_multi2 = GI.MultiPolygon([gi_poly2])

        
    # Convert to ArchGDAL geometries
    ag_poly1 = GI.convert(AG, gi_poly1)
    ag_poly2 = GI.convert(AG, gi_poly2)
    ag_multi1 = GI.convert(AG, gi_multi1)
    ag_multi2 = GI.convert(AG, gi_multi2)

    # Convert to LibGEOS geometries
    lg_poly1 = GI.convert(LG, gi_poly1)
    lg_poly2 = GI.convert(LG, gi_poly2)
    lg_multi1 = GI.convert(LG, gi_multi1)
    lg_multi2 = GI.convert(LG, gi_multi2)

    # Test matrix: all combinations (Polygons and MultiPolygons)
    geom_types = [
        ("GI.Polygon", gi_poly1, gi_poly2),
        ("GI.MultiPolygon", gi_multi1, gi_multi2),
        ("AG.Polygon", ag_poly1, ag_poly2),
        ("AG.MultiPolygon", ag_multi1, ag_multi2),
        ("LG.Polygon", lg_poly1, lg_poly2),
        ("LG.MultiPolygon", lg_multi1, lg_multi2),
    ]

    operations = [
        ("intersection", GO.intersection),
        ("union", GO.union),
        ("difference", GO.difference),
    ]

    for (op_name, op_func) in operations
        @testset "$op_name" begin
            for (type1_name, poly1a, poly1b) in geom_types
                for (type2_name, poly2a, poly2b) in geom_types
                    @testset "$type1_name x $type2_name" begin
                        result = _test_falseiferror(op_func, poly1a, poly2a; target = GI.PolygonTrait)
                        @test result
                    end
                end
            end
        end
    end

    @testset "Intersection with polygon with crs v/s without crs" begin
        @test GI.crs(first(GO.intersection(gi_poly_crs, gi_poly2; target = GI.PolygonTrait))) == GI.crs(gi_poly_crs)
        @test GI.crs(first(GO.intersection(ag_poly1, gi_poly_crs; target = GI.PolygonTrait))) == GI.crs(gi_poly_crs)
    end
end
