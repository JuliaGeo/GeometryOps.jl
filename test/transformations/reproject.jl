using Test
using GeoFormatTypes
import GeoInterface as GI
import GeometryOps as GO
import Proj
using ..TestHelpers

ring1 = GI.LinearRing([(1, 2), (7, 4), (5, 6), (1, 2)])
ring2 = GI.LinearRing([(11, 2), (20, 4), (15, 6), (11, 2)])
hole2 = GI.LinearRing([(14, 4), (16, 4), (17, 5), (14, 4)])

# Set up a regular tranformation of the points for reference
source_crs = convert(Proj.CRS, EPSG(4326))
target_crs = convert(Proj.CRS, EPSG(3857))
trans = Proj.Transformation(source_crs, target_crs; always_xy=true)

polygon1 = GI.Polygon([ring1])
polygon2 = GI.Polygon([ring2, hole2])
multipolygon = GI.MultiPolygon([polygon1, polygon2])

ref_points3857 = map(GI.getpoint(multipolygon)) do p
    trans([GI.x(p), GI.y(p)])
end

_xy(p) = GI.x(p), GI.y(p)

@testset_implementations "reproject" begin
    multipolygon3857 = GO.reproject($multipolygon, EPSG(4326), EPSG(3857))
    multipolygon4326 = GO.reproject($multipolygon; source_crs=EPSG(4326), target_crs=EPSG(4326))
    points4326_1 = collect(GI.getpoint($multipolygon))
    points4326_2 = collect.(GI.getcoord.(GI.getpoint(multipolygon4326)))
    points3857 = collect.(GI.getcoord.(GI.getpoint(multipolygon3857)))

    # Comparison to regular `trans` on points
    @test map(ref_points3857, points3857) do p1, p2
        all(map(isapprox, _xy(p1), _xy(p2))) 
    end |> all

    # Round trip comparison
    @test all(map((p1, p2) -> all(map(isapprox, _xy(p1), _xy(p2))), points4326_1, points4326_2))

    # Embedded crs check
    @test GI.crs(multipolygon3857) == EPSG(3857)
    @test GI.crs(multipolygon4326) == EPSG(4326)

    # Run it threaded over 100 replicates
    @test_nowarn GO.reproject([multipolygon3857 for _ in 1:100]; target_crs=EPSG(4326), threaded=true, calc_extent=true)

    utm32_wkt = """
    PROJCS["WGS 84 / UTM zone 32N",
        GEOGCS["WGS 84",
            DATUM["WGS_1984",
                SPHEROID["WGS 84",6378137,298.257223563,
                    AUTHORITY["EPSG","7030"]],
                AUTHORITY["EPSG","6326"]],
            PRIMEM["Greenwich",0,
                AUTHORITY["EPSG","8901"]],
            UNIT["degree",0.0174532925199433,
                AUTHORITY["EPSG","9122"]],
            AUTHORITY["EPSG","4326"]],
        PROJECTION["Transverse_Mercator"],
        PARAMETER["latitude_of_origin",0],
        PARAMETER["central_meridian",9],
        PARAMETER["scale_factor",0.9996],
        PARAMETER["false_easting",500000],
        PARAMETER["false_northing",0],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]],
        AXIS["Easting",EAST],
        AXIS["Northing",NORTH],
        AUTHORITY["EPSG","32632"]]
    """

    @test GO.reproject(multipolygon4326; source_crs="epsg:4326", target_crs="+proj=utm +zone=32 +datum=WGS84") ==
        GO.reproject(multipolygon4326; source_crs=EPSG(4326), target_crs=ProjString("+proj=utm +zone=32 +datum=WGS84 +type=crs")) ==
        GO.reproject(multipolygon4326; target_crs=EPSG(32632)) ==
        GO.reproject(multipolygon4326; target_crs="epsg:32632") ==
        GO.reproject(multipolygon4326; target_crs=utm32_wkt)

    GO.reproject(multipolygon4326; target_crs=ProjString("+proj=moll"))
end

