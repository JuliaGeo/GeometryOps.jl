using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import ArchGDAL as AG
import GeoFormatTypes as GFT
import GeometryBasics as GB
using GeometryOpsTestHelpers
using GeometryOps: Planar, Spherical, ChenMcMains, GEOS
using Random

#=
The native `buffer` against GEOS, whose buffer is the same JTS lineage — the
offset-curve primitives are transcribed from `OffsetSegmentGenerator`, so the
fillet vertices agree bit for bit and the two results differ only where GEOS
applies a heuristic we deliberately skip (see `src/methods/buffer_offset_curve.jl`).

The comparison is therefore a *relative symmetric-difference area*, at JTS
`BufferResultMatcher`'s own tolerance of 1e-3, with the observed agreement being
several orders tighter (`1e-17`-ish on these shapes). Vertex-level equality is
not the contract and is not asserted.
=#

const EMPTY_MP = LG.readgeom("MULTIPOLYGON EMPTY")

to_lg(g) = (GI.trait(g) isa GI.MultiPolygonTrait && GI.ngeom(g) == 0) ? EMPTY_MP :
           GI.convert(LG, g)

function rel_symdiff(actual, expected)
    a = to_lg(actual)
    ae = LG.area(expected)
    LG.isEmpty(a) && LG.isEmpty(expected) && return 0.0
    (LG.isEmpty(a) || LG.isEmpty(expected)) && return Inf
    return LG.area(LG.symmetricDifference(a, expected)) / ae
end

# The whole contract on one case: matches GEOS, and is a valid geometry.
function agrees_with_geos(geom, d; quadsegs = 8, atol = 1e-3)
    actual = GO.buffer(geom, d; quadsegs)
    expected = LG.buffer(GI.convert(LG, geom), d, quadsegs)
    lg = to_lg(actual)
    return LG.isValid(lg) && rel_symdiff(actual, expected) < atol
end

# ---------------------------------------------------------------------------

pt = GI.Point((0.3, -0.2))
mpt = GI.MultiPoint([(0.0, 0.0), (1.4, 0.0), (0.7, 1.2)])
line = GI.LineString([(0.0, 0.0), (3.0, 0.5), (4.0, 3.0), (1.5, 2.0)])
zigzag = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (4.0, 0.0)])
ring = GI.LinearRing([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)])
square = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
holed = GI.Polygon([
    [(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
    [(2.5, 0.8), (3.5, 0.8), (3.5, 3.2), (2.5, 3.2), (2.5, 0.8)],
])
concave = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (2.2, 1.0), (0.0, 4.0), (0.0, 0.0)]])
mpoly = GI.MultiPolygon([square, GI.Polygon([[(5.0, 0.0), (7.0, 0.0), (7.0, 2.0), (5.0, 0.0)]])])
collection = GI.GeometryCollection([pt, zigzag, concave])

@testset "Against LibGEOS" begin
    @testset "$name" for (name, geom, ds) in (
            ("point", pt, (1.0, 0.05)),
            ("multipoint", mpt, (1.0, 0.4)),
            ("line", line, (0.6, 2.0)),
            ("zigzag line", zigzag, (0.2, 1.5)),
            ("ring", ring, (0.5,)),
            ("square", square, (1.0, -0.5, -1.9)),
            ("polygon with a hole", holed, (0.3, 1.0, -0.3, -0.8, -1.5)),
            ("concave polygon", concave, (0.5, 2.0, -0.2)),
            ("multipolygon", mpoly, (0.3, 1.0, -0.3)),
            ("geometry collection", collection, (0.5, 2.0)))
        for d in ds
            @test agrees_with_geos(geom, d)
        end
    end

    @testset "quadsegs" begin
        for q in (1, 2, 4, 8, 16, 32)
            @test agrees_with_geos(pt, 1.0; quadsegs = q)
            @test agrees_with_geos(concave, 0.7; quadsegs = q)
            #-- the fillet segment count is JTS's, so the vertex count matches too
            @test GI.npoint(GO.buffer(pt, 1.0; quadsegs = q)) ==
                  GI.npoint(LG.buffer(GI.convert(LG, pt), 1.0, q))
        end
    end
end

#=
The same API over every geometry implementation GeometryOpsTestHelpers knows
about (GeoInterface wrappers always; ArchGDAL, GeometryBasics and LibGEOS via
their extensions), to check that nothing in the pipeline reads a concrete type.
=#
@testset_implementations "Implementations" begin
    @test GO.area(GO.buffer($square, 1.0)) ≈ GO.area(LG.buffer(GI.convert(LG, square), 1.0, 8)) rtol = 1e-9
    @test GI.trait(GO.buffer($square, -1.0)) isa GI.PolygonTrait
    @test GO.area(GO.buffer($line, 0.6)) ≈ GO.area(LG.buffer(GI.convert(LG, line), 0.6, 8)) rtol = 1e-9
end

@testset "Return shape" begin
    @test GI.trait(GO.buffer(pt, 1.0)) isa GI.PolygonTrait
    @test GI.trait(GO.buffer(line, 0.6)) isa GI.PolygonTrait
    #-- the most specific geometry: three overlapping discs merge into one polygon
    @test GI.trait(GO.buffer(mpt, 1.0)) isa GI.PolygonTrait
    #-- and stay separate when they do not
    @test GI.trait(GO.buffer(mpt, 0.2)) isa GI.MultiPolygonTrait
    @test GI.ngeom(GO.buffer(mpt, 0.2)) == 3
    #-- eroding the holed polygon splits it in two
    @test GI.trait(GO.buffer(holed, -0.5)) isa GI.MultiPolygonTrait
    #-- coordinates are always Float64 tuples, whatever went in
    @test first(GI.getpoint(GO.buffer(GI.Point(0, 0), 1))) isa Tuple{Float64, Float64}

    @testset "containers" begin
        #-- `apply` at `TraitTarget{GI.AbstractGeometryTrait}()`: containers come back
        v = GO.buffer([pt, line, square], 0.5)
        @test v isa Vector && length(v) == 3
        @test all(g -> GI.trait(g) isa GI.PolygonTrait, v)
        @test GI.trait(GO.buffer([pt pt; pt pt], 1.0)) === nothing  # a matrix, not a geometry
        fc = GI.FeatureCollection([GI.Feature(square; properties = (; a = 1)),
                                   GI.Feature(line; properties = (; a = 2))])
        bfc = GO.buffer(fc, 0.5)
        @test GI.trait(bfc) isa GI.FeatureCollectionTrait
        @test GI.nfeature(bfc) == 2
        @test GI.properties(GI.getfeature(bfc, 1)).a == 1
        @test GI.trait(GI.geometry(GI.getfeature(bfc, 2))) isa GI.PolygonTrait
    end

    @testset "crs and calc_extent" begin
        crs = GFT.EPSG(4326)
        withcrs = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]; crs)
        @test GI.crs(GO.buffer(withcrs, 0.2)) == crs
        @test GI.crs(GO.buffer(withcrs, -0.2)) == crs       # empty result keeps it too
        #-- `GI.extent` computes one on demand, so the EMBEDDED extent is the
        #-- thing `calc_extent` controls
        @test getfield(GO.buffer(square, 1.0), :extent) === nothing
        ext = getfield(GO.buffer(square, 1.0; calc_extent = true), :extent)
        @test ext isa GI.Extents.Extent
        @test all(ext.X .≈ (-1.0, 5.0)) && all(ext.Y .≈ (-1.0, 5.0))
    end
end

@testset "Degenerate cases" begin
    @testset "empty results" begin
        #-- negative distance on anything of dimension < 2
        for g in (pt, mpt, line, zigzag, ring, GI.GeometryCollection([pt, line]))
            for d in (-1.0, -1e-9)
                b = GO.buffer(g, d)
                @test GI.trait(b) isa GI.MultiPolygonTrait
                @test GI.ngeom(b) == 0
                @test GI.npoint(b) == 0
            end
        end
        #-- total erosion of an areal input. This is the `buffer(linestring, -1.0)`
        #-- `BoundsError` from the old GEOS-forwarding default: `tuples` -> `rebuild`
        #-- called `first(child_geoms)` on an empty vector.
        @test GI.ngeom(GO.buffer(square, -3.0)) == 0
        @test GI.ngeom(GO.buffer(holed, -3.0)) == 0
        @test GO.area(GO.buffer(square, -2.0)) == 0.0   # exactly the inradius
    end

    @testset "zero distance" begin
        #-- areal input at d == 0 is the input; everything else is empty, as in GEOS
        @test GO.area(GO.buffer(square, 0.0)) ≈ 16.0
        @test GO.area(GO.buffer(holed, 0.0)) ≈ GO.area(holed)
        for g in (pt, mpt, line, ring)
            @test GI.ngeom(GO.buffer(g, 0.0)) == 0
        end
    end

    @testset "degenerate linework" begin
        #-- a "ring" with fewer than 3 distinct vertices buffers as a line,
        #-- a line with one distinct point buffers as a point (JTS `getRingCurve`)
        collapsed = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)]])
        @test agrees_with_geos(collapsed, 0.5)
        repeated = GI.LineString([(0.0, 0.0), (0.0, 0.0), (1.0, 0.0), (1.0, 0.0)])
        @test agrees_with_geos(repeated, 0.4)
        allsame = GI.LineString([(2.0, 2.0), (2.0, 2.0), (2.0, 2.0)])
        @test GO.area(GO.buffer(allsame, 0.5)) ≈ GO.area(GO.buffer(GI.Point(2.0, 2.0), 0.5))
        #-- an exact 180° reversal, which needs the half-disc fillet
        hairpin = GI.LineString([(0.0, 0.0), (3.0, 0.0), (0.0, 0.0)])
        @test agrees_with_geos(hairpin, 0.5)
        #-- a closed line is linework, not an area: it buffers on both sides
        closed = GI.LineString([(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)])
        @test agrees_with_geos(closed, 0.5)
        #-- non-finite coordinates are rejected, not propagated
        @test_throws ArgumentError GO.buffer(GI.Point(NaN, 0.0), 1.0)
    end

    #=
    JTS `BufferCurveSetBuilder.isRingFullyEroded`: a ring whose own erosion is
    already empty emits no curve. The answers below are what the engine returns
    without the pre-check too — it is a speed guard, not a correctness one — so
    each case also asserts the answer, not just that something was dropped.
    =#
    @testset "fully eroded rings" begin
        function ngon(n, R = 10.0)
            v = [(R * cos(2π * i / n), R * sin(2π * i / n)) for i in 0:(n - 1)]
            return GI.Polygon([push!(v, v[1])])
        end
        #-- shell: an n-gon eroded past its apothem is empty, and cheap. Without the
        #-- pre-check this same call builds ~n^2/2 self-crossings to reach `∅`.
        big = ngon(400)                              # apothem 9.99969, envelope 20 x 20
        @test GI.ngeom(GO.buffer(big, -10.5)) == 0
        @test GO._raw_offset_curves(Tuple{Float64, Float64}, big, -10.5, 8) == []
        #-- ... and the near miss, past collapse but inside the envelope test, still
        #-- goes through the engine and still comes back empty
        @test 2 * 9.9999 < 20.0                      # envelope is 2R = 20, so no drop
        @test !isempty(GO._raw_offset_curves(Tuple{Float64, Float64}, big, -9.9999, 8))
        @test GI.ngeom(GO.buffer(big, -9.9999)) == 0
        #-- ... and just short of collapse the result is a real, non-empty region
        @test GO.area(GO.buffer(big, -9.9)) > 0
        @test agrees_with_geos(ngon(40), -9.0)
        #-- the test is `2|d| >= envMin`, where JTS has `>`: at equality the erosion
        #-- is a measure-zero set and the engine returns nothing for it anyway
        @test GO._raw_offset_curves(Tuple{Float64, Float64}, square, -2.0, 8) == []
        @test GO.area(GO.buffer(square, -2.0)) == 0.0

        #-- a triangle uses the exact incircle radius, not the envelope
        tri = GI.Polygon([[(0.0, 0.0), (12.0, 0.0), (0.0, 9.0), (0.0, 0.0)]])
        inr = 2 * 54 / (12 + 9 + 15)                 # 2·area / perimeter = 3
        @test GI.ngeom(GO.buffer(tri, -1.001 * inr)) == 0
        @test GO.area(GO.buffer(tri, -0.5 * inr)) ≈ 54 / 4 rtol = 1e-9
        #-- the envelope test would not have fired here: min(12, 9) > 2·3.003
        @test 2 * 1.001 * inr < 9.0

        #-- a dropped shell takes its holes with it: the hole's curve alone would
        #-- have nothing to be a hole in
        holed_thin = GI.Polygon([
            [(0.0, 0.0), (20.0, 0.0), (20.0, 1.0), (0.0, 1.0), (0.0, 0.0)],
            [(5.0, 0.4), (5.0, 0.6), (15.0, 0.6), (15.0, 0.4), (5.0, 0.4)],
        ])
        @test GI.ngeom(GO.buffer(holed_thin, -0.6)) == 0

        #-- a hole the dilation swallows is dropped, and the shell fills it
        narrow_hole = GI.Polygon([
            [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0), (0.0, 0.0)],
            [(5.0, 9.5), (5.0, 10.5), (15.0, 10.5), (15.0, 9.5), (5.0, 9.5)],
        ])
        @test agrees_with_geos(narrow_hole, 0.51)
        @test GI.nhole(GO.buffer(narrow_hole, 0.51)) == 0
        #-- the near miss: 2d is exactly the hole's width, so the hole is kept
        @test agrees_with_geos(narrow_hole, 0.5)
        @test agrees_with_geos(narrow_hole, 0.3)
        @test GI.nhole(GO.buffer(narrow_hole, 0.3)) == 1

        #-- a MultiPolygon part is dropped on its own, the others survive
        mp = GI.MultiPolygon([
            GI.Polygon([[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]]),
            GI.Polygon([[(20.0, 0.0), (30.0, 0.0), (30.0, 1.0), (20.0, 1.0), (20.0, 0.0)]]),
        ])
        @test agrees_with_geos(mp, -0.6)
        @test GO.area(GO.buffer(mp, -0.6)) ≈ 8.8^2 rtol = 1e-9
    end

    #=
    `prune_eroded_rings = false` takes the pre-check out and reads the engine's
    own answer at collapse. The two must agree everywhere — the check is a speed
    guard over an identity that is already unconditional — so this is the test
    that keeps the keyword honest on both sides of every threshold it owns.
    =#
    @testset "the fully-eroded pre-check is optional" begin
        ngon(n, R = 10.0) = GI.Polygon([push!(
            [(R * cos(2π * i / n), R * sin(2π * i / n)) for i in 0:(n - 1)],
            (R * 1.0, 0.0))])
        shell = ngon(60)                             # apothem 9.9863, envelope 20 x 20
        #-- across collapse: inside it, at it, and well past it
        for d in (-9.9, -9.986, -9.99, -10.0, -10.5, -12.0)
            on = GO.buffer(shell, d)
            off = GO.buffer(shell, d; prune_eroded_rings = false)
            @test GI.ngeom(to_lg(on)) == GI.ngeom(to_lg(off))
            @test GO.area(on) ≈ GO.area(off) atol = 1e-12
        end
        #-- the hole side of the check: `d > 0` swallowing a 1-wide hole
        narrow = GI.Polygon([
            [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0), (0.0, 0.0)],
            [(5.0, 9.5), (5.0, 10.5), (15.0, 10.5), (15.0, 9.5), (5.0, 9.5)],
        ])
        for d in (0.3, 0.5, 0.51, 0.6)
            @test GO.area(GO.buffer(narrow, d)) ≈
                  GO.area(GO.buffer(narrow, d; prune_eroded_rings = false)) atol = 1e-9
        end
        #-- degenerate linework is dropped whatever the keyword says: it has no
        #-- interior to erode, and the generator has no curve for it
        collapsed = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)]])
        @test GI.ngeom(GO.buffer(collapsed, -0.5; prune_eroded_rings = false)) == 0
        #-- and it is the curve set, not the answer, that the keyword changes
        T = Tuple{Float64, Float64}
        @test GO._raw_offset_curves(T, shell, -10.5, 8) == []
        @test !isempty(GO._raw_offset_curves(T, shell, -10.5, 8;
                                             prune_eroded_rings = false))
    end

    #=
    P10. `_offset_ring_orientation` decides which side of a ring the curve goes
    on, so a wrong answer buffers the complement rather than perturbing the
    result. The shoelace it replaced summed products of absolute coordinates: a
    clockwise unit square at (1e8, 1e8) has terms of size 1e16 and an area of 1,
    so the sum came out exactly zero, the shell was never reversed, and the
    buffer had area 0.04 against GEOS's 1.43.
    =#
    @testset "translated coordinates" begin
        unitsq(b) = [(b, b), (b + 1.0, b), (b + 1.0, b + 1.0), (b, b + 1.0)]
        #-- `GO.area` is itself an absolute-coordinate shoelace, so the result is
        #-- brought back to the origin before it is measured
        shift(g, dx) = GO.apply(GI.PointTrait(), g) do p
            (GI.x(p) + dx, GI.y(p) + dx)
        end
        at_origin = GO.area(GO.buffer(GI.Polygon([[unitsq(0.0)..., (0.0, 0.0)]]), 0.1))
        for base in (1e8, 1e10), rev in (false, true)
            v = rev ? reverse(unitsq(base)) : unitsq(base)
            g = GI.Polygon([[v..., v[1]]])
            @test agrees_with_geos(g, 0.1)
            #-- translation invariance, to the grid the far coordinates are on
            @test GO.area(shift(GO.buffer(g, 0.1), -base)) ≈ at_origin rtol = 1e-4
        end
        #-- a hole at the same magnitude, whose orientation is decided separately
        holed_far = GI.Polygon([
            [(1e8, 1e8), (1e8 + 10, 1e8), (1e8 + 10, 1e8 + 10), (1e8, 1e8 + 10), (1e8, 1e8)],
            [(1e8 + 3, 1e8 + 3), (1e8 + 3, 1e8 + 7), (1e8 + 7, 1e8 + 7), (1e8 + 7, 1e8 + 3),
             (1e8 + 3, 1e8 + 3)],
        ])
        @test agrees_with_geos(holed_far, 0.5)
        @test agrees_with_geos(holed_far, -0.5)
    end

    #=
    P15. A non-finite distance used to reach the generator and surface as an
    `InexactError` from a NaN; a distance below the coordinate grid used to
    return the input by accident, depending on where the rounding landed.
    =#
    @testset "distance guard" begin
        for d in (Inf, -Inf, NaN)
            @test_throws ArgumentError GO.buffer(square, d)
            @test_throws ArgumentError GO.buffer(line, d)
        end
        base = 1e10
        v = [(base, base), (base + 1, base), (base + 1, base + 1), (base, base + 1)]
        far = GI.Polygon([[v..., v[1]]])
        grid = 0.5 * eps(base + 1)                   # half an ulp there, ~9.5e-7
        back(g) = GO.apply(GI.PointTrait(), g) do p  # `GO.area` at 1e10 is useless
            (GI.x(p) - base, GI.y(p) - base)
        end
        #-- below the grid: no offset point can leave its own grid point, so the
        #-- distance is zero and an areal geometry comes back unchanged
        @test GO.area(back(GO.buffer(far, 0.9 * grid))) == 1.0
        @test GO.area(back(GO.buffer(far, -0.9 * grid))) == 1.0
        #-- above it, a real buffer
        @test GO.area(back(GO.buffer(far, 1e-3))) > 1.003
        #-- lower dimensions are empty there, exactly as at `d == 0`
        @test GI.ngeom(GO.buffer(line, 0.4 * eps(4.0))) == 0
        @test GI.ngeom(GO.buffer(pt, 0.4 * eps(0.3))) == 0
    end

    @testset "self-intersecting input" begin
        #-- invalid input has no defined answer, but must not throw or produce
        #-- invalid output; GEOS agrees on this one because both follow the
        #-- curve's local turn sense
        bowtie = GI.Polygon([[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0), (0.0, 0.0)]])
        @test LG.isValid(to_lg(GO.buffer(bowtie, 0.3)))
        crossing = GI.LineString([(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)])
        @test agrees_with_geos(crossing, 0.4)
    end
end

@testset "Algorithm and keywords" begin
    @test GO.buffer(ChenMcMains(), square, 1.0) == GO.buffer(square, 1.0)
    @test GO.buffer(Planar(), square, 1.0) == GO.buffer(square, 1.0)
    @test GO.manifold(ChenMcMains()) === Planar()
    @test GO.rebuild(ChenMcMains(; quadsegs = 3), Planar()).quadsegs == 3
    @test ChenMcMains(; prune_eroded_rings = false).prune_eroded_rings == false
    @test GO.rebuild(ChenMcMains(; prune_eroded_rings = false),
                     Planar()).prune_eroded_rings == false

    @testset "unsupported styles error, naming GEOS" begin
        for kw in ((; endCapStyle = :flat), (; endCapStyle = :square),
                   (; joinStyle = :mitre), (; joinStyle = :bevel))
            e = try GO.buffer(square, 1.0; kw...) catch e; e end
            @test e isa ArgumentError
            @test occursin("GEOS", e.msg)
        end
        #-- the round defaults are accepted, so a call that spelled them out survives
        @test GO.buffer(square, 1.0; endCapStyle = :round, joinStyle = :round,
                        mitreLimit = 2.0, single_sided = false) == GO.buffer(square, 1.0)
        @test_throws ArgumentError ChenMcMains(; quadsegs = 0)
        #-- `single_sided` is rejected too, and its message names the value
        e = try GO.buffer(square, 1.0; single_sided = true) catch e; e end
        @test e isa ArgumentError
        @test occursin("single_sided = true", e.msg) && occursin("GEOS", e.msg)
    end

    @testset "spherical is not implemented" begin
        @test_throws ArgumentError ChenMcMains(Spherical())
        @test_throws ArgumentError GO.buffer(Spherical(), square, 1.0)
    end

    @testset "the GEOS path still works" begin
        @test GO.equals(GO.buffer(GEOS(), square, 1.0), LG.buffer(GI.convert(LG, square), 1.0))
        @test rel_symdiff(GO.buffer(GEOS(; joinStyle = :mitre), square, 1.0),
                          LG.bufferWithStyle(GI.convert(LG, square), 1.0;
                                             joinStyle = LG.GEOSBUF_JOIN_MITRE)) < 1e-12
        @test_nowarn GO.buffer(GEOS(), [square mpoly], 1.0)
        @test_nowarn GO.buffer(GEOS(), collection, 1.0)
    end
end

#=
A small differential fuzz: seeded star polygons (jittered even angles, so always
valid) and random-walk polylines, at log-uniform distances of both signs. Sized
to stay inside this file's ~30 s budget; the standing 536-case corpus and the
5000-case raw-curve fuzz live outside the test suite.
=#
@testset "Fuzz vs LibGEOS" begin
    rng = Xoshiro(20260917)
    nfail = 0
    for i in 1:80
        geom, d = if isodd(i)
            n = rand(rng, 4:14)
            base = range(0, 2π; length = n + 1)[1:n]
            θ = base .+ rand(rng, n) .* (2π / n) .* 0.7
            r = 1 .+ 4 .* rand(rng, n)
            pts = [(r[k] * cos(θ[k]), r[k] * sin(θ[k])) for k in 1:n]
            (GI.Polygon([[pts..., pts[1]]]), (rand(rng) < 0.5 ? -1 : 1) * exp(1.2 * randn(rng) - 1.5))
        else
            n = rand(rng, 2:12)
            x, y, a = 0.0, 0.0, rand(rng) * 2π
            pts = [(0.0, 0.0)]
            for _ in 2:n
                a += rand(rng) < 0.3 ? 3 * (2 * rand(rng) - 1) : 0.8 * randn(rng)
                x += cos(a); y += sin(a)
                push!(pts, (x, y))
            end
            (GI.LineString(pts), exp(1.2 * randn(rng) - 0.5))
        end
        #-- skip inputs GEOS itself calls invalid: there is no reference answer
        GI.trait(geom) isa GI.PolygonTrait && !LG.isValid(GI.convert(LG, geom)) && continue
        ok = try
            agrees_with_geos(geom, d)
        catch e
            @error "buffer threw" i geom d exception = e
            false
        end
        ok || (nfail += 1)
    end
    @test nfail == 0
end
