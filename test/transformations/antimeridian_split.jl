using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import ArchGDAL as AG
using GeometryOpsTestHelpers

const SPH = GO.Spherical()

# ---------------------------------------------------------------------------
# helpers

# every near-pole vertex run (|lat| ≈ 90) in a piece, in ring order
function pole_rows(piece)
    rows = Vector{Vector{Tuple{Float64,Float64}}}()
    for r in GI.getring(piece)
        pts = [(Float64(GI.x(q)), Float64(GI.y(q))) for q in GI.getpoint(r)]
        i = 1
        while i <= length(pts)
            if abs(pts[i][2]) >= 89.999999
                j = i
                while j < length(pts) && abs(pts[j+1][2]) >= 89.999999
                    j += 1
                end
                push!(rows, pts[i:j])
                i = j + 1
            else
                i += 1
            end
        end
    end
    return rows
end

npole_verts(piece) = sum(r -> count(q -> abs(Float64(GI.y(q))) >= 89.999999, GI.getpoint(r)),
                         GI.getring(piece); init = 0)

# global gates on an emitted MultiPolygon: high-level LibGEOS validity, no
# |Δlon| > 180 anywhere, all longitudes in the closed branch [λn-360, λn], and
# (optionally) spherical area conservation vs a reference.
function assert_gates(mp; λn = 180.0, check_jump = true, area_ref = nothing, rtol = 1e-12)
    @test GI.trait(mp) isa GI.MultiPolygonTrait
    for p in GI.getgeom(mp)
        @test LG.isValid(GI.convert(LG, p))
        for r in GI.getring(p)
            lons = [Float64(GI.x(q)) for q in GI.getpoint(r)]
            @test all(l -> λn - 360.0 - 1e-9 <= l <= λn + 1e-9, lons)
            if check_jump
                for i in 1:(length(lons)-1)
                    @test abs(lons[i+1] - lons[i]) <= 180.0 + 1e-9
                end
            end
        end
    end
    if area_ref !== nothing
        a_out = sum(p -> GO.area(SPH, p), GI.getgeom(mp); init = 0.0)
        @test isapprox(a_out, area_ref; rtol)
    end
    return nothing
end

# orient a single-ring polygon to the smaller of its two spherical regions (the
# bounded interpretation, e.g. the polar cap rather than its global complement)
function orient_small(pts)
    a = GI.Polygon([pts]); b = GI.Polygon([reverse(pts)])
    GO.area(SPH, a) <= GO.area(SPH, b) ? a : b
end

# a lat-constant 12-gon, oriented to the *small* cap (encloses the near pole)
capring(lat) = (r = [(Float64(l), lat) for l in 15.0:30.0:345.0]; push!(r, r[1]); r)
mkcap(lat) = orient_small(capring(lat))

piece_lons(p) = sort(unique(round.([Float64(GI.x(q)) for q in GI.getpoint(GI.getexterior(p))]; digits = 6)))

# ---------------------------------------------------------------------------
@testset "simple box crossing the seam -> 2 branch pieces" begin
    box = GI.Polygon([[(170.0, 40.0), (-170.0, 40.0), (-170.0, 50.0), (170.0, 50.0), (170.0, 40.0)]])
    mp = GO.antimeridian_split(box)
    @test GI.ngeom(mp) == 2
    assert_gates(mp; area_ref = GO.area(SPH, box))

    lonsets = sort([piece_lons(p) for p in GI.getgeom(mp)])
    @test lonsets == sort([[170.0, 180.0], [-180.0, -170.0]])   # west lip 180, east lip -180
    # each piece lives strictly within one 360-wide branch (max lon span < 180)
    for p in GI.getgeom(mp)
        lons = [Float64(GI.x(q)) for q in GI.getpoint(GI.getexterior(p))]
        @test maximum(lons) - minimum(lons) <= 180.0
    end
end

# ---------------------------------------------------------------------------
@testset "south/north caps enclosing a pole -> 1 piece with pole row" begin
    for (lat, polelat, lip_south) in ((-70.0, -90.0, true), (70.0, 90.0, false))
        cap = mkcap(lat)
        mp = GO.antimeridian_split(cap; pole_spacing = 5.0)
        @test GI.ngeom(mp) == 1
        assert_gates(mp; area_ref = GO.area(SPH, cap))

        rows = pole_rows(GI.getgeom(mp, 1))
        @test length(rows) == 1
        row = rows[1]
        @test length(row) == 73                         # ns=72 -> 71 interior + 2 corners
        lons = first.(row)
        d = diff(lons)
        @test all(<(0), d) || all(>(0), d)              # strictly monotone
        @test maximum(abs.(d)) <= 5.0 + 1e-9            # step <= spacing
        # both corners exactly on the seam lips at the pole
        @test (180.0, polelat) in row
        @test (-180.0, polelat) in row
    end
end

# ---------------------------------------------------------------------------
@testset "pole_spacing = nothing -> corners only, area identical" begin
    cap = mkcap(-70.0)
    mp_on = GO.antimeridian_split(cap; pole_spacing = 5.0)
    mp_off = GO.antimeridian_split(cap; pole_spacing = nothing)

    row_off = pole_rows(GI.getgeom(mp_off, 1))[1]
    @test length(row_off) == 2                          # just the two corners
    @test (180.0, -90.0) in row_off && (-180.0, -90.0) in row_off

    # the infill points are the same kernel point (zero-area) -> areas are equal
    a_on = sum(p -> GO.area(SPH, p), GI.getgeom(mp_on))
    a_off = sum(p -> GO.area(SPH, p), GI.getgeom(mp_off))
    @test a_on == a_off
end

# ---------------------------------------------------------------------------
@testset "wedge with a vertex AT the pole (not enclosing) -> no pole row" begin
    wedge = orient_small([(100.0, -80.0), (140.0, -80.0), (120.0, -90.0), (100.0, -80.0)])
    mp = GO.antimeridian_split(wedge)
    @test GI.ngeom(mp) == 1
    assert_gates(mp; area_ref = GO.area(SPH, wedge))
    p = GI.getgeom(mp, 1)
    @test all(r -> length(r) == 1, pole_rows(p))        # a single pole vertex, not a resampled row
    @test npole_verts(p) == 1
end

# ---------------------------------------------------------------------------
@testset "touch-without-crossing and no-crossing fast path" begin
    touch = GI.Polygon([[(170.0, 0.0), (180.0, 10.0), (170.0, 20.0), (160.0, 10.0), (170.0, 0.0)]])
    mp = GO.antimeridian_split(touch)
    @test GI.ngeom(mp) == 1
    assert_gates(mp; area_ref = GO.area(SPH, touch))

    # a box around lon 0 cannot interact with the ±180 seam -> fast path, verbatim
    b0 = GI.Polygon([[(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0), (-10.0, -10.0)]])
    mp = GO.antimeridian_split(b0)
    @test GI.ngeom(mp) == 1
    got = [(Float64(GI.x(q)), Float64(GI.y(q))) for q in GI.getpoint(GI.getexterior(GI.getgeom(mp, 1)))]
    want = [(Float64(GI.x(q)), Float64(GI.y(q))) for q in GI.getpoint(GI.getexterior(b0))]
    @test got == want                                   # output ≡ input
end

# ---------------------------------------------------------------------------
@testset "holes" begin
    # hole crossing the seam (through the same arrangement machinery)
    outer = [(160.0, -20.0), (160.0, 20.0), (-160.0, 20.0), (-160.0, -20.0), (160.0, -20.0)]
    hole = [(175.0, -5.0), (175.0, 5.0), (-175.0, 5.0), (-175.0, -5.0), (175.0, -5.0)]
    holed = GI.Polygon([outer, hole])
    mp = GO.antimeridian_split(holed)
    @test GI.ngeom(mp) == 2
    assert_gates(mp; area_ref = GO.area(SPH, holed))

    # hole entirely on one side of the seam
    outer2 = [(150.0, -20.0), (150.0, 20.0), (-150.0, 20.0), (-150.0, -20.0), (150.0, -20.0)]
    hole2 = [(160.0, -5.0), (160.0, 5.0), (170.0, 5.0), (170.0, -5.0), (160.0, -5.0)]   # west of seam
    holed2 = GI.Polygon([outer2, hole2])
    mp2 = GO.antimeridian_split(holed2)
    @test GI.ngeom(mp2) == 2
    assert_gates(mp2; area_ref = GO.area(SPH, holed2))
end

# ---------------------------------------------------------------------------
@testset "arbitrary seam longitude" begin
    box = GI.Polygon([[(170.0, 40.0), (-170.0, 40.0), (-170.0, 50.0), (170.0, 50.0), (170.0, 40.0)]])

    # default == explicit 180.0, bit-for-bit
    m_def = GO.antimeridian_split(box)
    m_180 = GO.antimeridian_split(box; antimeridian = 180.0)
    @test collect(GI.getpoint(m_def)) == collect(GI.getpoint(m_180))

    # a Greenwich-crossing box splits at antimeridian = 0 (lips 0 / -360), NOT at default
    gbox = GI.Polygon([[(-10.0, 40.0), (10.0, 40.0), (10.0, 50.0), (-10.0, 50.0), (-10.0, 40.0)]])
    m0 = GO.antimeridian_split(gbox; antimeridian = 0.0)
    @test GI.ngeom(m0) == 2
    assert_gates(m0; λn = 0.0, area_ref = GO.area(SPH, gbox))
    lonsets = sort([piece_lons(p) for p in GI.getgeom(m0)])
    @test lonsets == sort([[-10.0, 0.0], [-360.0, -350.0]])   # west lip 0, east lip -360
    @test GI.ngeom(GO.antimeridian_split(gbox)) == 1          # untouched at default seam

    # Pacific seam (330° == -30° normalised): a box crossing lon -30 splits
    pbox = GI.Polygon([[(-40.0, 10.0), (-20.0, 10.0), (-20.0, 20.0), (-40.0, 20.0), (-40.0, 10.0)]])
    mpac = GO.antimeridian_split(pbox; antimeridian = 330.0)
    @test GI.ngeom(mpac) == 2
    assert_gates(mpac; λn = -30.0, area_ref = GO.area(SPH, pbox))
    # -30 and 330 normalise to the same seam
    @test collect(GI.getpoint(GO.antimeridian_split(pbox; antimeridian = -30.0))) ==
          collect(GI.getpoint(mpac))
end

# ---------------------------------------------------------------------------
@testset "rotated pole" begin
    # a cap around geographic (0,0); with the rotated pole placed at (0,0) the
    # cap encloses the rotated north pole and crosses the rotated seam.
    capg = Tuple{Float64,Float64}[(20.0 * cosd(az), 20.0 * sind(az)) for az in 0.0:30.0:330.0]
    push!(capg, capg[1])
    cap = GI.Polygon([capg])

    mp = GO.antimeridian_split(cap; north_pole = (0.0, 0.0))
    @test GI.ngeom(mp) == 1
    # output is in the ROTATED frame (default seam), so gates hold there; area is
    # rotation-invariant, so conservation is against the input area directly.
    assert_gates(mp; area_ref = GO.area(SPH, cap))
    rows = pole_rows(GI.getgeom(mp, 1))
    @test length(rows) == 1                             # pole row at the rotated ±90
    @test length(rows[1]) >= 3

    # rotation is an isometry -> area conserved to near machine precision
    a_out = sum(p -> GO.area(SPH, p), GI.getgeom(mp))
    @test isapprox(a_out, GO.area(SPH, cap); rtol = 1e-12)
end

# ---------------------------------------------------------------------------
@testset "unsupported traits throw" begin
    @test_throws ArgumentError GO.antimeridian_split(GI.LineString([(0.0, 0.0), (1.0, 1.0)]))
    @test_throws ArgumentError GO.antimeridian_split(GI.Point(0.0, 0.0))
end

# ---------------------------------------------------------------------------
# Ported literal fixtures from the Python `antimeridian` package's test suite
# (inlined — no network). These overlap the synthetic cases but pin external
# parity for the two canonical shapes.
@testset "python antimeridian package fixtures" begin
    # `crossing` fixture: a quad straddling the seam -> 2 pieces
    crossing = GI.Polygon([[(170.0, -10.0), (-170.0, -10.0), (-170.0, 10.0), (170.0, 10.0), (170.0, -10.0)]])
    mc = GO.antimeridian_split(crossing)
    @test GI.ngeom(mc) == 2
    assert_gates(mc; area_ref = GO.area(SPH, crossing))

    # `north-pole` fixture: a band enclosing the north pole -> 1 piece + pole row
    npoly = mkcap(85.0)
    mn = GO.antimeridian_split(npoly)
    @test GI.ngeom(mn) == 1
    assert_gates(mn; area_ref = GO.area(SPH, npoly))
    @test length(pole_rows(GI.getgeom(mn, 1))) == 1
end

# ---------------------------------------------------------------------------
# implementation-generic: the split works on any GeoInterface polygon
impl_box = GI.Polygon([[(170.0, 40.0), (-170.0, 40.0), (-170.0, 50.0), (170.0, 50.0), (170.0, 40.0)]])
@testset_implementations "crossing box across implementations" begin
    mp = GO.antimeridian_split($impl_box)
    @test GI.trait(mp) isa GI.MultiPolygonTrait
    @test GI.ngeom(mp) == 2
    @test all(p -> LG.isValid(GI.convert(LG, p)), GI.getgeom(mp))
    a_out = sum(p -> GO.area(SPH, p), GI.getgeom(mp))
    @test isapprox(a_out, GO.area(SPH, $impl_box); rtol = 1e-12)
end

# ---------------------------------------------------------------------------
# Natural Earth real data (gated on data availability, same pattern as the other
# overlayng NE tests — NaturalEarth is a test dependency).
ne_ok = false
ne_names = String[]; ne_geoms = Any[]
try
    import NaturalEarth, GeoJSON
    fc = NaturalEarth.naturalearth("admin_0_countries", 110)
    for f in fc
        g = GeoJSON.geometry(f)
        (g === nothing || GI.npoint(g) == 0) && continue
        nm = try; string(f.NAME); catch; "?"; end
        push!(ne_names, nm); push!(ne_geoms, GO.tuples(g))
    end
    global ne_ok = length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
end

@testset "Natural Earth antimeridian split" begin
    if !ne_ok
        @test_skip "Natural Earth data unavailable"
    else
        byname(nm) = ne_geoms[findfirst(==(nm), ne_names)]
        for (nm, expect_pieces, wants_polerow) in
                (("Russia", 14, false), ("Fiji", 3, false), ("Antarctica", 8, true))
            g = byname(nm)
            mp = GO.antimeridian_split(g)
            @test GI.ngeom(mp) == expect_pieces
            assert_gates(mp; area_ref = GO.area(SPH, g), rtol = 1e-13)
            haspolerow = any(p -> !isempty(pole_rows(p)), GI.getgeom(mp))
            @test haspolerow == wants_polerow
        end

        # rotated-pole NE sanity: split Russia about a Pacific pole. Output is in
        # the rotated frame; validity and (isometry-preserved) area still hold.
        # (Antarctica is deliberately NOT used here: NE encodes its polar cap with
        # an antimeridian slit that is only OGC-valid because the slit lies on the
        # ±180 seam — rotating it turns the slit into an interior degeneracy, so a
        # rotated-frame split of Antarctica is legitimately invalid.)
        ru = byname("Russia")
        mpr = GO.antimeridian_split(ru; north_pole = (-150.0, 10.0))
        assert_gates(mpr; area_ref = GO.area(SPH, ru), rtol = 1e-12)
    end
end
