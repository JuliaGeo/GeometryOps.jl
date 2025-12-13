using Test

import GeometryOps as GO, GeoInterface as GI
import TGGeometry

# This test file is only really certifying that the extension is working...
# we rely on TGGeometry.jl's tests to verify correctness, since it's not an exact
# library and a lot of our tests do check exactness.

# That's why it isn't included in the main polygon test suite.

point = (0.0, 0.0)
multipoint = GI.MultiPoint([(0.0, 0.0), (1.0, 1.0)])
linestring = GI.LineString([(0.0, 0.0), (1.0, 1.0)])
multilinestring = GI.MultiLineString([[(0.0, 0.0), (1.0, 1.0)], [(2.0, 2.0), (3.0, 3.0)]])
polygon = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
multipolygon = GI.MultiPolygon([polygon, GO.transform(p -> (GI.x(p), GI.y(p) + 1.0), polygon)])

_xplus5(g) = GO.transform(p -> (GI.x(p) + 5.0, GI.y(p)), g)

disjoint_point = _xplus5(point)
disjoint_multipoint = _xplus5(multipoint)
disjoint_linestring = _xplus5(linestring)
disjoint_multilinestring = _xplus5(multilinestring)
disjoint_polygon = _xplus5(polygon)
disjoint_multipolygon = _xplus5(multipolygon)

@testset "Internal consistency with TGGeometry.jl" begin
    for funsym in TGGeometry.TG_PREDICATES
        for geom1 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
            for geom2 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
                @testset let predicate = funsym, geom1 = geom1, geom2 = geom2
                    @test Base.getproperty(TGGeometry, predicate)(geom1, geom2) == Base.getproperty(GO, predicate)(GO.TG(), geom1, geom2)
                end
            end
            for geom2 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
                @testset let predicate = funsym, geom1 = geom1, geom2 = geom2
                    @test Base.getproperty(TGGeometry, predicate)(geom1, geom2) == Base.getproperty(GO, predicate)(GO.TG(), geom1, geom2)
                end
            end
        end
    end
end

@testset "Consistency with GeometryOps algorithms for simple cases" begin
    for funsym in TGGeometry.TG_PREDICATES
        for geom1 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
            for geom2 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
                @testset let predicate = funsym, geom1 = geom1, geom2 = geom2
                    @test Base.getproperty(TGGeometry, predicate)(geom1, geom2) == Base.getproperty(GO, predicate)(GO.TG(), geom1, geom2)
                end
            end
            for geom2 in (point, multipoint, linestring, multilinestring, polygon, multipolygon)
                @testset let predicate = funsym, geom1 = geom1, geom2 = geom2
                    @test Base.getproperty(TGGeometry, predicate)(geom1, geom2) == Base.getproperty(GO, predicate)(GO.TG(), geom1, geom2)
                end
            end
        end
    end
end