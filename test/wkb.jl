using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import Random
import NaturalEarth
using GeometryOpsTestHelpers
using GeometryOpsTestHelpers: parse_wkb, write_wkb

# LibGEOS is the reference codec.  High-level API only (raw GEOS C calls abort the process).
const CTX = LG.get_global_context()
lg_wkb(g) = LG.writegeom(GI.convert(LG, g), LG.WKBWriter(CTX), CTX)

# Count vertices in a (multi)polygon, for ns/vertex reporting.
function nverts(geom)
    t = GI.trait(geom)
    if t isa GI.PolygonTrait
        return sum(GI.npoint, GI.getgeom(geom); init = 0)
    elseif t isa GI.MultiPolygonTrait
        return sum(nverts, GI.getgeom(geom); init = 0)
    else
        return 0
    end
end

# ------------------------------------------------------------------
# Structure-aware little-endian -> big-endian converter, used to build
# big-endian and mixed-endianness fixtures from our own writer output.
# (Test-only helper; correctness over speed.)
mutable struct _Cur
    buf::Vector{UInt8}
    pos::Int   # 1-based
end
_r_u8!(c) = (v = c.buf[c.pos]; c.pos += 1; v)
_r_u32le!(c) = (v = reinterpret(UInt32, c.buf[c.pos:c.pos+3])[1]; c.pos += 4; v)
_r_8!(c) = (v = c.buf[c.pos:c.pos+7]; c.pos += 8; v)
_be32(v::UInt32) = reverse(reinterpret(UInt8, [v]))

function _swap_geom!(c::_Cur, out::Vector{UInt8})
    order = _r_u8!(c)
    @assert order == 0x01 "converter expects little-endian input"
    push!(out, 0x00)                       # emit big-endian
    typ = _r_u32le!(c)
    append!(out, _be32(typ))
    code = typ % 1000
    if code == 3
        nrings = _r_u32le!(c); append!(out, _be32(nrings))
        for _ in 1:nrings
            npts = _r_u32le!(c); append!(out, _be32(npts))
            for _ in 1:(npts * 2)
                append!(out, reverse(_r_8!(c)))
            end
        end
    elseif code == 6
        ngeoms = _r_u32le!(c); append!(out, _be32(ngeoms))
        for _ in 1:ngeoms
            _swap_geom!(c, out)
        end
    else
        error("test converter only handles polygon/multipolygon")
    end
    return out
end
to_be(le::Vector{UInt8}) = _swap_geom!(_Cur(le, 1), UInt8[])

# ------------------------------------------------------------------
# Fixtures
unit_square = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]])
two_holes = GI.Polygon([
    [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
    [(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0), (1.0, 1.0)],
    [(5.0, 5.0), (5.0, 6.0), (6.0, 6.0), (6.0, 5.0), (5.0, 5.0)],
])
mp_holes = GI.MultiPolygon([
    [[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)],
     [(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0), (2.0, 2.0)]],
    [[(20.0, 20.0), (23.0, 20.0), (23.0, 23.0), (20.0, 20.0)]],
])
# Coordinates that stress floating-point round-tripping.
stress_vals = [
    (1 / 3, 1e-300),
    (-0.0, 0.0),
    (nextfloat(0.0), floatmax(Float64)),
    (prevfloat(0.0), 1.7976931348623157e308),
    (1 / 3, 1e-300),
]
stress_poly = GI.Polygon([stress_vals])

@testset "Round trip: fixtures (both directions)" begin
    for (name, g) in [
        ("unit square", unit_square),
        ("two holes", two_holes),
        ("multipolygon with holes", mp_holes),
        ("float stress", stress_poly),
    ]
        @testset "$name" begin
            # LibGEOS-written WKB -> parse_wkb == GO.tuples(g)
            @test parse_wkb(lg_wkb(g)) == GO.tuples(g)
            # write_wkb(g) -> LibGEOS parse, coordinate-exact
            @test GO.tuples(LG.readgeom(write_wkb(g))) == GO.tuples(g)
            # our writer round trip
            @test parse_wkb(write_wkb(g)) == GO.tuples(g)
            # our writer is byte-identical to LibGEOS' (LE ISO)
            @test write_wkb(g) == lg_wkb(g)
        end
    end
end

@testset "Bit-level float exactness" begin
    g = parse_wkb(write_wkb(stress_poly))
    pts = collect(GI.getpoint(GI.getexterior(g)))
    for (pt, v) in zip(pts, stress_vals)
        @test reinterpret(UInt64, GI.x(pt)) == reinterpret(UInt64, v[1])
        @test reinterpret(UInt64, GI.y(pt)) == reinterpret(UInt64, v[2])
    end
    # And through LibGEOS (binary, so exact)
    g2 = parse_wkb(lg_wkb(stress_poly))
    pts2 = collect(GI.getpoint(GI.getexterior(g2)))
    for (pt, v) in zip(pts2, stress_vals)
        @test reinterpret(UInt64, GI.x(pt)) == reinterpret(UInt64, v[1])
        @test reinterpret(UInt64, GI.y(pt)) == reinterpret(UInt64, v[2])
    end
end

@testset "Empty geometries" begin
    # Empty multipolygon (LibGEOS can represent this one).
    emp_wkb = lg_wkb(GI.MultiPolygon([[[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]]]))
    emp = parse_wkb(UInt8[1, 6, 0, 0, 0, 0, 0, 0, 0])
    @test GI.trait(emp) isa GI.MultiPolygonTrait
    @test GI.ngeom(emp) == 0
    @test write_wkb(emp) == UInt8[1, 6, 0, 0, 0, 0, 0, 0, 0]
    @test parse_wkb(write_wkb(emp)) == emp
    # LibGEOS agrees an empty multipolygon is 9 bytes.
    empty_mp_lg = LG.readgeom(UInt8[1, 6, 0, 0, 0, 0, 0, 0, 0])
    @test LG.writegeom(empty_mp_lg, LG.WKBWriter(CTX), CTX) == UInt8[1, 6, 0, 0, 0, 0, 0, 0, 0]

    # Empty polygon (LibGEOS' GI.convert can't build this, so test our codec directly).
    ep = parse_wkb(UInt8[1, 3, 0, 0, 0, 0, 0, 0, 0])
    @test GI.trait(ep) isa GI.PolygonTrait
    @test GI.ngeom(ep) == 0
    @test write_wkb(ep) == UInt8[1, 3, 0, 0, 0, 0, 0, 0, 0]
    @test parse_wkb(write_wkb(ep)) == ep
end

@testset "Degenerate but parseable rings" begin
    # Rings with too few points to be valid, and unclosed — still valid WKB.
    for coords in ([(1.0, 2.0)], [(1.0, 2.0), (3.0, 4.0)], Tuple{Float64,Float64}[])
        rings = Vector{Tuple{Float64,Float64}}[coords]
        w = UInt8[]
        push!(w, 0x01)
        append!(w, reinterpret(UInt8, [htol(UInt32(3))]))
        append!(w, reinterpret(UInt8, [htol(UInt32(1))]))          # 1 ring
        append!(w, reinterpret(UInt8, [htol(UInt32(length(coords)))]))
        for (x, y) in coords
            append!(w, reinterpret(UInt8, [htol(reinterpret(UInt64, x))]))
            append!(w, reinterpret(UInt8, [htol(reinterpret(UInt64, y))]))
        end
        g = parse_wkb(w)
        @test GI.npoint(GI.getexterior(g)) == length(coords)
        @test write_wkb(g) == w    # exact re-serialization
    end
end

@testset "Big-endian and mixed endianness" begin
    for g in (unit_square, two_holes, mp_holes)
        le = write_wkb(g)
        be = to_be(le)
        @test be[1] == 0x00                     # confirm we built a big-endian buffer
        @test parse_wkb(be) == parse_wkb(le)
        @test parse_wkb(be) == GO.tuples(g)
    end

    # Mixed-endianness multipolygon: sub-polygon 1 little-endian, sub-polygon 2 big-endian.
    poly1 = GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]])
    poly2 = GI.Polygon([[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 5.0)]])
    mixed = UInt8[]
    push!(mixed, 0x01)                                        # LE multipolygon header
    append!(mixed, reinterpret(UInt8, [htol(UInt32(6))]))
    append!(mixed, reinterpret(UInt8, [htol(UInt32(2))]))
    append!(mixed, write_wkb(poly1))                         # LE sub-polygon
    append!(mixed, to_be(write_wkb(poly2)))                  # BE sub-polygon
    expected = GI.MultiPolygon([
        [[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]],
        [[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 5.0)]],
    ])
    @test parse_wkb(mixed) == GO.tuples(expected)
end

@testset "EWKB SRID accepted and ignored" begin
    g = GI.Polygon([[(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (1.0, 2.0)]])
    le = write_wkb(g)
    typ = reinterpret(UInt32, le[2:5])[1]
    out = UInt8[le[1]]
    append!(out, reinterpret(UInt8, [htol(typ | 0x20000000)]))   # set SRID flag
    append!(out, reinterpret(UInt8, [htol(UInt32(4326))]))       # SRID payload
    append!(out, le[6:end])
    @test parse_wkb(out) == GO.tuples(g)
end

@testset "Z / M dimensions rejected" begin
    g = GI.Polygon([[(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (1.0, 2.0)]])
    le = write_wkb(g)
    # ISO PolygonZ (1003)
    iso = copy(le); iso[2:5] = reinterpret(UInt8, [htol(UInt32(1003))])
    @test_throws ArgumentError parse_wkb(iso)
    # ISO PolygonM (2003) and PolygonZM (3003)
    for code in (2003, 3003)
        b = copy(le); b[2:5] = reinterpret(UInt8, [htol(UInt32(code))])
        @test_throws ArgumentError parse_wkb(b)
    end
    # EWKB Z / M flag bits
    for flag in (0x80000000, 0x40000000, 0xc0000000)
        b = copy(le); b[2:5] = reinterpret(UInt8, [htol(UInt32(3) | flag)])
        @test_throws ArgumentError parse_wkb(b)
    end
end

@testset "Unsupported types and malformed buffers throw cleanly" begin
    # Point / LineString / MultiPoint / MultiLineString / GeometryCollection
    for code in (1, 2, 4, 5, 7)
        b = UInt8[0x01]
        append!(b, reinterpret(UInt8, [htol(UInt32(code))]))
        append!(b, reinterpret(UInt8, [htol(UInt32(0))]))
        @test_throws ArgumentError parse_wkb(b)
    end
    # empty buffer / truncated header / truncated count / truncated coords
    @test_throws ArgumentError parse_wkb(UInt8[])
    @test_throws ArgumentError parse_wkb(UInt8[1, 3])
    @test_throws ArgumentError parse_wkb(UInt8[1, 3, 0, 0, 0])          # no ring count
    @test_throws ArgumentError parse_wkb(UInt8[1, 3, 0, 0, 0, 1, 0, 0, 0])  # ring count but no point count
    # absurd counts must not OOM / segfault
    @test_throws ArgumentError parse_wkb(UInt8[1, 3, 0, 0, 0, 0xff, 0xff, 0xff, 0xff])
    @test_throws ArgumentError parse_wkb(UInt8[1, 3, 0, 0, 0, 1, 0, 0, 0, 0xff, 0xff, 0xff, 0xff])
    @test_throws ArgumentError parse_wkb(UInt8[1, 6, 0, 0, 0, 0xff, 0xff, 0xff, 0xff])
    # bad byte-order flag
    @test_throws ArgumentError parse_wkb(UInt8[2, 3, 0, 0, 0, 0, 0, 0, 0])
    # multipolygon whose element is not a polygon
    bad_mp = UInt8[0x01]
    append!(bad_mp, reinterpret(UInt8, [htol(UInt32(6))]))
    append!(bad_mp, reinterpret(UInt8, [htol(UInt32(1))]))
    append!(bad_mp, UInt8[0x01])
    append!(bad_mp, reinterpret(UInt8, [htol(UInt32(1))]))   # a Point inside the MultiPolygon
    append!(bad_mp, reinterpret(UInt8, [htol(reinterpret(UInt64, 0.0))]))
    append!(bad_mp, reinterpret(UInt8, [htol(reinterpret(UInt64, 0.0))]))
    @test_throws ArgumentError parse_wkb(bad_mp)
    # write_wkb rejects non-polygonal geometry
    @test_throws ArgumentError write_wkb(GI.Point(1.0, 2.0))
    @test_throws ArgumentError write_wkb(GI.LineString([(0.0, 0.0), (1.0, 1.0)]))
end

@testset "Seeded fuzz vs LibGEOS" begin
    Random.seed!(0xC0FFEE)
    rndcoord() = rand() * 200 - 100
    function rndring()
        k = rand(3:7)
        pts = [(rndcoord(), rndcoord()) for _ in 1:k]
        push!(pts, pts[1])   # close the ring, keep it valid for GEOS (>= 4 pts)
        return pts
    end
    rndpoly() = GI.Polygon([rndring() for _ in 1:rand(1:3)])
    rndmpoly() = GI.MultiPolygon([[rndring() for _ in 1:rand(1:3)] for _ in 1:rand(1:3)])

    for _ in 1:250
        g = rndpoly()
        @test parse_wkb(write_wkb(g)) == g                         # our codec self-consistency
        @test parse_wkb(lg_wkb(g)) == GO.tuples(g)                 # LibGEOS write -> our parse
        @test GO.tuples(LG.readgeom(write_wkb(g))) == GO.tuples(g) # our write -> LibGEOS parse
    end
    for _ in 1:150
        g = rndmpoly()
        @test parse_wkb(write_wkb(g)) == g
        @test parse_wkb(lg_wkb(g)) == GO.tuples(g)
        @test GO.tuples(LG.readgeom(write_wkb(g))) == GO.tuples(g)
    end
end

@testset "Natural Earth round trip + benchmark" begin
    countries = try
        NaturalEarth.naturalearth("admin_0_countries", 110)
    catch e
        @info "Skipping Natural Earth WKB tests (data unavailable)" exception = e
        nothing
    end
    if countries !== nothing
        geoms = [g for g in countries.geometry if GI.trait(g) isa Union{GI.PolygonTrait,GI.MultiPolygonTrait}]
        @test length(geoms) > 100

        lgbufs = [lg_wkb(g) for g in geoms]
        for (g, buf) in zip(geoms, lgbufs)
            @test parse_wkb(buf) == GO.tuples(g)                        # both directions,
            @test GO.tuples(LG.readgeom(write_wkb(g))) == GO.tuples(g)  # coordinate-exact
        end

        # --- Benchmark: ns/vertex, ours vs LibGEOS.readgeom + GO.tuples ---
        totalverts = sum(nverts, geoms)
        tupgeoms = [parse_wkb(b) for b in lgbufs]   # tuple geometries, the harness' write input

        parse_mine() = (for b in lgbufs; parse_wkb(b); end; nothing)
        parse_lg() = (for b in lgbufs; GO.tuples(LG.readgeom(b)); end; nothing)
        wr = LG.WKBWriter(CTX)
        lggeoms = [GI.convert(LG, g) for g in geoms]
        write_mine() = (for g in tupgeoms; write_wkb(g); end; nothing)
        write_lg() = (for lg in lggeoms; LG.writegeom(lg, wr, CTX); end; nothing)

        parse_mine(); parse_lg(); write_mine(); write_lg()   # warmup / compile
        tpm = minimum(@elapsed(parse_mine()) for _ in 1:10)
        tpl = minimum(@elapsed(parse_lg()) for _ in 1:10)
        twm = minimum(@elapsed(write_mine()) for _ in 1:10)
        twl = minimum(@elapsed(write_lg()) for _ in 1:10)

        @info "WKB benchmark (Natural Earth 110m countries)" nfeatures = length(geoms) totalverts
        @info "parse ns/vertex" parse_wkb = tpm / totalverts * 1e9 libgeos_readgeom_tuples = tpl / totalverts * 1e9 speedup = tpl / tpm
        @info "write ns/vertex" write_wkb = twm / totalverts * 1e9 libgeos_wkbwriter = twl / totalverts * 1e9 speedup = twl / twm
        @test tpm > 0   # sanity: benchmark actually ran
    end
end
