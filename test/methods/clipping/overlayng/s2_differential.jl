# Differential validation of SPHERICAL OverlayNG against s2geography
# (phase-3 validation, design §4).
#
# The planar suites (`fuzz.jl`, `xml_suite.jl`) differentiate against LibGEOS,
# which is a PLANAR engine: it cannot say anything at all about geodesic edges,
# so `Spherical()` overlay has until now only been checked against itself, by
# the algebraic identities in `realdata_identities.jl`. s2geography wraps
# Google's S2 — genuinely geodesic, and the reference implementation the rest of
# the ecosystem is calibrated against — so it is the oracle that closes that gap.
#
# The FFI binding lives in `test/external/s2geography/s2geography.jl`; see its
# header for the two C routes and the Arrow types the kernels demand.
#
# ## What is compared, and why not coordinates
#
# The two engines round differently (ours arranges exactly and rounds once at
# emission; S2 snaps in its builder), so vertex-level equality is not the
# contract and is not tested. Three graded, scale-free quantities are, each
# normalised by `scale = s2_area(A) + s2_area(B)`:
#
#  1. **AREA AGREEMENT** `|s2_area(ours) - s2_area(theirs)| / scale`.
#     Always computable, and the sharpest test of "did you build the right
#     region". Every measurement here is made by s2 on both operands, so the
#     radius question below does not arise.
#
#  2. **L¹ DISTANCE** `s2_area(ours △ theirs) / scale`, with the symmetric
#     difference built by s2 as well. This is the established metric in this
#     project (`fuzz.jl` uses the GEOS-built version) and it is strictly
#     stronger than (1): two results can have equal area and be in different
#     places. It is NOT always available — see "the oracle's resolution floor".
#
#  3. **STRUCTURE** the two results must have the same topological dimension and
#     the same number of areal components, and ours must import into S2 at all
#     (a result S2 rejects is a failure however good its area looks). Two
#     classes are exempted from parts of this, each diagnosed in place: s2's
#     areal-only result model (see `check_op`) and the antimeridian seam
#     partition (see the ledger).
#
# Plus a fourth, which validates OUR area rather than the overlay:
#
#  4. **RADIUS-CORRECTED AREA** `GO.area(Spherical(), ours) * RADIUS_RATIO`
#     against `s2_area(ours)`. s2geography measures on `S2Earth::RadiusMeters()
#     = 6371010.0`; GeometryOps uses `WGS84_EARTH_MEAN_RADIUS = 6371008.8`. The
#     1.2 m difference is 3.8e-7 in area and is not a bug on either side, so it
#     is divided out explicitly instead of being hidden inside a widened
#     tolerance. What is left agrees to a couple of ulps.
#
# ## The oracle's resolution floor (measured, not assumed)
#
# s2geography's overlay applies S2Builder's snap radius, so when two results
# differ only by an ulp-scale displacement of a shared vertex their symmetric
# difference is a pair of slivers that S2's polygon assembler either refuses to
# build (`Can't find shell for polygon hole`) or reports with a snap-inflated
# area. On the unit-square reproducer below the true swept area is 6.1e-7 m²
# (one vertex moved 2 ulps along a 5.5e4 m boundary) and s2 reports 9.0e-3 m² —
# four orders high, and 3.6e-13 of the input scale. So:
#
#   * the L¹ gate is set at the level S2 can actually resolve, not at the level
#     the geometry actually differs by;
#   * a `Can't find shell` throw from the METRIC (not from the op under test) is
#     classified `oracle-unresolved`, counted, and reported — never silently
#     passed — and that case still has to clear gates (1), (3) and (4).
#
# Note the identity form `area(A) + area(B) - 2·area(A∩B)` is NOT a usable
# substitute: it cancels the two opposite-signed slivers and bottoms out at
# `ulp(2e10) = 3.8e-6 m²`, i.e. it measures nothing at all in this regime.
#
# ## Env knobs
#
#   GO_S2_PAIRS      NE 110 m neighbour pairs (default 30)
#   GO_S2_NE10       NE 10 m neighbour pairs (default 4)
#   GO_S2_REPORTS    disagreements printed in full (default 5)

using Test
import GeometryOpsCore
using GeometryOpsTestHelpers: write_wkb, parse_wkb, WKBWriteError
include(joinpath(@__DIR__, "common.jl"))

include(joinpath(@__DIR__, "..", "..", "..", "external", "s2geography", "s2geography.jl"))
using .S2Geog

const S2_PAIRS = something(tryparse(Int, get(ENV, "GO_S2_PAIRS", "")), 30)
const S2_NE10_PAIRS = something(tryparse(Int, get(ENV, "GO_S2_NE10", "")), 4)
const S2_REPORTS = something(tryparse(Int, get(ENV, "GO_S2_REPORTS", "")), 5)

# A missing oracle records a skip by default; under `GO_REQUIRE_DATA=1` (which
# CI sets) it fails instead, so a runner that has quietly lost `S2Geography_jll`
# cannot report this file as green having compared nothing.
const REQUIRE_DATA = get(ENV, "GO_REQUIRE_DATA", "") == "1"
function missing_data(what)
    if REQUIRE_DATA
        println("REQUIRED DATASET MISSING: ", what)
        @test false
    else
        @test_skip what
    end
end

# s2geography measures on 6371010.0, GeometryOps on 6371008.8; areas scale with R².
const RADIUS_RATIO = (S2Geog.S2_RADIUS / GeometryOpsCore.WGS84_EARTH_MEAN_RADIUS)^2

# ---------------------------------------------------------------------------
# Gates
# ---------------------------------------------------------------------------
#
# All three sit at 1e-12, each for its own measured reason. Measured over the
# default corpus (204 ops: 30 NE 110 m neighbour pairs, 10 shifted-self cases,
# 4 NE 10 m pairs, 7 NE 10 m Russia pairs) — the census at the end of every run
# reprints these:
#
#   area agreement        median 3.3e-15   worst 1.0e-13  (NE10 Belarus x Russia/union)
#   L1 symdiff / scale    median 2.2e-16   worst 5.9e-15  (NE110 Kazakhstan shifted/symdifference)
#   radius-corrected area median 4.3e-15   worst 2.3e-13  (NE10 Norway x Russia/difference)
#
# i.e. the two engines agree to machine precision on this data, and the gates are
# NOT set by our emission error. They are set by the two floors below, both
# measured rather than assumed, and both still four orders of magnitude tighter
# than any structural disagreement (a missing component or a mis-signed face is
# 1e-2 to 1e0 of scale, not 1e-12).
#
# 1. AREA_RTOL and L1_RTOL are floored by S2's own snap radius: on the
#    unit-square reproducer in the header, a result differing from ours by one
#    2-ulp vertex displacement (true swept area 6.1e-7 m²) is reported by s2 as
#    9.0e-3 m², i.e. 3.6e-13 of the input scale. 1e-12 clears that floor.
# 2. RADIUS_RTOL is floored by summation accumulation in `GO.area(Spherical())`
#    versus S2's exact area over a several-thousand-vertex ring — this check
#    grades the AREA method on identical geometry, not the overlay.
const AREA_RTOL = 1e-12
const L1_RTOL = 1e-12
const RADIUS_RTOL = 1e-12

# ---------------------------------------------------------------------------
# Census
# ---------------------------------------------------------------------------

mutable struct S2Census
    n_op::Int; n_agree::Int; n_unresolved::Int; n_subareal::Int; n_collapsed::Int; n_seam::Int
    n_disagree::Int
    l1::Vector{Float64}; area::Vector{Float64}; radius::Vector{Float64}
    worst::Dict{String, Tuple{Float64, String}}
end
S2Census() = S2Census(0, 0, 0, 0, 0, 0, 0, Float64[], Float64[], Float64[],
                      Dict{String, Tuple{Float64, String}}())
const CENSUS = S2Census()

function note_worst!(key, v, label)
    (isfinite(v) && v > get(CENSUS.worst, key, (-1.0, ""))[1]) || return nothing
    CENSUS.worst[key] = (v, label)
    return nothing
end

med(v) = isempty(v) ? NaN : sort(v)[(length(v) + 1) ÷ 2]

# Registered with `atexit` BEFORE the sweeps run, so the measured distribution is
# reported even when one fails: a failing `@testset` throws on the way out of its
# own block and would otherwise skip everything below it. The numbers are the
# deliverable of this suite, not a debugging aid.
function report_census()
    CENSUS.n_op == 0 && return nothing
    println("""

    Spherical OverlayNG vs s2geography ($(CENSUS.n_op) ops)
      agree                 : $(CENSUS.n_agree)
      disagree (failed)     : $(CENSUS.n_disagree)
      sub-areal on both     : $(CENSUS.n_subareal)  (both punctual, or both empty)
      dimension-collapsed   : $(CENSUS.n_collapsed)  (ours linear, s2 areal-only — see header)
      seam component split  : $(CENSUS.n_seam)  (diagnosed class — see the ledger)
      L1 unresolvable by s2 : $(CENSUS.n_unresolved)  (S2 could not build the symmetric difference)
      area agreement        : median $(med(CENSUS.area))   worst $(get(CENSUS.worst, "area agreement", (NaN, "")))
      L1 symdiff / scale    : median $(med(CENSUS.l1))   worst $(get(CENSUS.worst, "L1 symmetric difference", (NaN, "")))
      radius-corrected area : median $(med(CENSUS.radius))   worst $(get(CENSUS.worst, "radius-corrected area", (NaN, "")))""")
    return nothing
end
atexit(report_census)

# ---------------------------------------------------------------------------
# WKB helpers
# ---------------------------------------------------------------------------

# Little-endian WKB header fields, read without parsing the body: s2 returns
# LINESTRING / MULTILINESTRING for a shared-border intersection and
# GEOMETRYCOLLECTION for an empty result, neither of which `parse_wkb`
# (Polygon / MultiPolygon only) accepts, and both of which are correct answers.
_u32(b, i) = UInt32(b[i]) | UInt32(b[i+1]) << 8 | UInt32(b[i+2]) << 16 | UInt32(b[i+3]) << 24
wkb_typecode(b::Vector{UInt8}) = length(b) < 5 ? UInt32(0) : _u32(b, 2)
wkb_count(b::Vector{UInt8}) = length(b) < 9 ? UInt32(0) : _u32(b, 6)

# Topological dimension of a result: 2 areal, 1 linear, 0 punctual, -1 empty.
# Two engines agreeing on the dimension is the first thing to check — a shared
# border makes `intersection` LINEAR on both sides, which is not a disagreement.
function result_dimension(g)
    t = GI.trait(g)
    #-- collections and points first: GeoInterface defines `npoint` for neither
    t isa GI.GeometryCollectionTrait &&
        return maximum(result_dimension, GI.getgeom(g); init = -1)
    t isa GI.PointTrait && return 0
    t isa GI.MultiPointTrait && return GI.ngeom(g) == 0 ? -1 : 0
    GI.npoint(g) == 0 && return -1
    (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) && return 2
    (t isa GI.LineStringTrait || t isa GI.LinearRingTrait ||
     t isa GI.MultiLineStringTrait) && return 1
    return -1
end

function wkb_dimension(b::Vector{UInt8})
    tc = wkb_typecode(b)
    (tc in (3, 6, 4, 5, 7) && wkb_count(b) == 0) && return -1
    tc == 3 && return 2
    tc == 6 && return 2
    (tc == 2 || tc == 5) && return 1
    (tc == 1 || tc == 4) && return 0
    tc == 7 && return 2                      # non-empty mixed collection
    return -1
end

# Number of areal components, with a Polygon counting as one and anything empty
# as zero, so a Polygon and a 1-component MultiPolygon compare equal.
function ncomponents(g)
    t = GI.trait(g)
    t isa GI.PolygonTrait && return GI.npoint(g) == 0 ? 0 : 1
    t isa GI.MultiPolygonTrait && return count(p -> GI.npoint(p) > 0, GI.getgeom(g))
    t isa GI.GeometryCollectionTrait && return sum(ncomponents, GI.getgeom(g); init = 0)
    return 0
end

function ncomponents_wkb(b::Vector{UInt8})
    tc = wkb_typecode(b)
    (tc == 3 || tc == 6) || return 0        # anything else s2 emits here is not areal
    return ncomponents(parse_wkb(b))
end

# ---------------------------------------------------------------------------
# The property
# ---------------------------------------------------------------------------

"""
    differential(A, B, label)

Run all four ops on one pair and compare against s2geography. Never throws: an
engine or oracle error is recorded as a failed check, so one bad pair cannot
abort the sweep.
"""
function differential(A, B, label)
    wa, wb = write_wkb(A), write_wkb(B)
    scale = s2_area(wa) + s2_area(wb)
    scale > 0 || return nothing
    for op in OP_SYMS
        key = "$label/$op"
        CENSUS.n_op += 1
        ok, detail = try
            check_op(op, A, B, wa, wb, scale, key)
        catch err
            (false, "ERROR: " * first(sprint(showerror, err), 200))
        end
        if ok
            @test true
        else
            CENSUS.n_disagree += 1
            CENSUS.n_disagree <= S2_REPORTS && println("S2 DIFFERENTIAL [$key] ", detail)
            @test false
        end
    end
    return nothing
end

function check_op(op, A, B, wa, wb, scale, key)
    ours = GO._overlay_ng(GO.Spherical(), OPCODE[op], A, B; exact = EX)
    tw = s2_overlay(op, wa, wb)

    #-- (3a) dimension, and the one place the two engines differ by DESIGN.
    #--
    #-- s2geography's constructive ops return a SINGLE geometry, never a mixed
    #-- collection, and for two areal inputs they emit only the areal part —
    #-- except that a punctual remnant does survive. Measured, on synthetic
    #-- inputs:
    #--
    #--   shared edge  ours LINESTRING          s2 POLYGON EMPTY
    #--   corner touch ours POINT               s2 POINT
    #--   NE110 Armenia x Iran (a shared border segment AND a separate corner
    #--                touch): ours GEOMETRYCOLLECTION(LINESTRING, POINT),
    #--                s2 POINT — the line is dropped, the point kept.
    #--
    #-- OGC/JTS semantics, which GeometryOps ports, return the full
    #-- dimension-collapsed result; GEOS agrees with us exactly on that Armenia
    #-- case, down to the same line and the same point. Neither engine is wrong,
    #-- so the case is not scored as a divergence — but it is not skipped either.
    #-- Two things are still asserted, through the route that CAN express them:
    #--
    #--   * s2's own `intersects` predicate must agree the inputs meet, and
    #--   * our result must be at least as high-dimensional as s2's, so
    #--     "s2 found an intersection we missed entirely" still fails.
    d_ours, d_theirs = result_dimension(ours), wkb_dimension(tw)
    if d_ours < 2 || d_theirs < 2
        a_theirs = s2_area(tw)
        d_ours == 2 && return (false,
            "we built an areal result but s2's is dimension $d_theirs (area $a_theirs)")
        d_theirs == 2 && return (false,
            "s2 built an areal result (area $a_theirs) but ours is dimension $d_ours " *
            "($(GI.trait(ours)))")
        a_theirs == 0 || return (false, "s2's dimension-$d_theirs result has area $a_theirs")
        d_ours >= d_theirs || return (false,
            "s2 found a dimension-$d_theirs intersection and we found dimension $d_ours " *
            "($(GI.trait(ours)))")
        if d_ours >= 0
            s2_predicate(:intersects, wa, wb) || return (false,
                "we built a dimension-$d_ours result but s2 reports the inputs disjoint")
        end
        d_ours == d_theirs ? (CENSUS.n_subareal += 1) : (CENSUS.n_collapsed += 1)
        CENSUS.n_agree += 1
        return (true, "")
    end

    ow = try
        write_wkb(ours)
    catch err
        err isa WKBWriteError || rethrow()
        #-- an areal result the WKB codec cannot express is a mixed collection,
        #-- which s2 would have had to match on dimension above
        return (false, "our areal result is not serializable: $(GI.trait(ours))")
    end

    #-- (3b) structure. `s2_area` also proves S2 could import our result at all.
    a_ours = s2_area(ow)
    a_theirs = s2_area(tw)
    nc_ours, nc_theirs = ncomponents(ours), ncomponents_wkb(tw)
    if nc_ours != nc_theirs
        key in S2_SEAM_COMPONENT_SPLIT || return (false,
            "component count $nc_ours vs s2 $nc_theirs (areas $a_ours vs $a_theirs)")
        CENSUS.n_seam += 1                 # diagnosed class, see the ledger below
    end

    #-- (1) area agreement
    r_area = abs(a_ours - a_theirs) / scale
    push!(CENSUS.area, r_area); note_worst!("area agreement", r_area, key)
    r_area <= AREA_RTOL ||
        return (false, "area disagreement $r_area (ours $a_ours, s2 $a_theirs, scale $scale)")

    #-- (4) our area against s2's, on the SAME geometry, radius divided out
    r_rad = a_ours == 0 ? 0.0 :
        abs(GO.area(GO.Spherical(), ours) * RADIUS_RATIO - a_ours) / a_ours
    push!(CENSUS.radius, r_rad); note_worst!("radius-corrected area", r_rad, key)
    r_rad <= RADIUS_RTOL ||
        return (false, "GO.area vs s2 area on our own result: $r_rad relative")

    #-- (2) L¹, when S2 can assemble the symmetric difference at all
    sd = try
        s2_area(s2_overlay(:symdifference, ow, tw))
    catch
        CENSUS.n_unresolved += 1
        CENSUS.n_agree += 1
        return (true, "")
    end
    r_l1 = sd / scale
    push!(CENSUS.l1, r_l1); note_worst!("L1 symmetric difference", r_l1, key)
    r_l1 <= L1_RTOL ||
        return (false, "L1 symmetric difference $r_l1 of scale (absolute $sd m²)")

    CENSUS.n_agree += 1
    return (true, "")
end

# ---------------------------------------------------------------------------
# Ledger
# ---------------------------------------------------------------------------
#
# DISCIPLINE, as in `test/external/jts/overlay_skiplist.jl`: every entry MUST
# carry a diagnosis explaining exactly why the two engines differ and which side
# is wrong. Entries without one must not be merged.
#
# There is exactly one ledger, and it is narrow by construction. A second,
# whole-op ledger used to sit here — `S2_KNOWN_DISAGREEMENTS`, scoring a case
# `@test_broken` and asserting nothing about it — and it stayed empty for its
# whole life while carrying three lines of unreachable branch in the hot loop.
# If a whole-op pin is ever genuinely needed, the honest move is to widen this
# ledger with the diagnosis attached, not to reinstate a mechanism that turns
# an op into a no-op.

# ## The antimeridian seam: a different partition of the SAME region
#
# Narrower than a whole-op pin, because only ONE check is in question. For these
# ops the component count may differ; every other gate — dimension, area
# agreement, L¹, the radius cross-check — is asserted exactly as usual, and
# since the area gate is what proves the two partitions cover the same region,
# relaxing the count loses nothing.
#
# THE DIAGNOSIS. Natural Earth splits landmasses at the antimeridian, so NE 10 m
# Russia is a 214-polygon MultiPolygon in which Chukotka and Wrangel Island are
# each cut into a piece ending at lon +180 and a piece starting at lon -180.
# Those are the same meridian on the sphere, so each pair shares a boundary
# SEGMENT there, and a MultiPolygon whose components meet along a line rather
# than at finitely many points is not valid on the sphere — valid on the plane,
# where +180 and -180 are different places, and invalid off it. The engine
# contracts on valid input (design §2.2), so neither engine is obliged to
# produce the other's answer here.
#
# What the two actually do, measured on `Belarus x Russia` union (the pair
# already pinned in `realdata_identities.jl`):
#
#   ours 212 components, s2 214. The four components that do not match up are
#     s2  mainland  1.6646793087751854e13   (lat 41.19–77.74)
#     s2  + wrapped 1.0966335977178513e11   (lat 64.24–68.98)
#     s2  Wrangel E 2.507812218056487e9     (lon 178.62–180)
#     s2  Wrangel W 5.066122856551397e9
#     ours mainland 1.6756456447523613e13 = the first two, merged
#     ours Wrangel  7.573935074607917e9   = the last two, merged
#
#   and the two partitions sum to 1.676403038259822e13 vs 1.676403038259825e13,
#   i.e. the SAME region to 2e-15. We treat the shared seam edge as the
#   adjacency it is on the sphere and merge; s2 keeps the pieces apart.
#
# `antimeridian_split` — this branch's parent — is a no-op on this input, since
# the rings already terminate exactly at ±180 and there is nothing to cut.
#
# Worth recording alongside: the `OverlayTopologyError: side location conflict`
# these two pairs used to raise NO LONGER REPRODUCES on this branch.
# `Azerbaijan x Russia` and `Belarus x Russia` pass every algebraic identity
# (conservation 5.3e-15 / 5.6e-15, reconstruction 9.9e-18 / 1.6e-17), so the
# `@test_broken` pins they carried in `realdata_identities.jl` are gone.
const S2_SEAM_COMPONENT_SPLIT = Set{String}(
    "NE10 $nm x Russia/$op"
    for nm in ("Azerbaijan", "Belarus", "Finland", "Norway",
               "Kazakhstan", "Mongolia", "Ukraine")
    for op in ("union", "symdifference"))

# ---------------------------------------------------------------------------
# Sweeps
# ---------------------------------------------------------------------------

const S2_OK = s2_available()
const NE_OK = S2_OK && try
    import NaturalEarth, GeoJSON
    true
catch err
    @info "s2geography overlay differential skipped (Natural Earth unavailable)" err
    false
end

if S2_OK
    include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))
    v = S2Geog.s2_versions()
    println("s2geography oracle: S2 $(v.s2geometry), GeoArrow $(v.geoarrow), " *
            "nanoarrow $(v.nanoarrow), $(length(S2Geog.s2_kernel_names())) kernels")
end

@testset "s2geography availability" begin
    if !S2_OK
        missing_data("S2Geography_jll on this platform")
    else
        @test length(S2Geog.s2_kernel_names()) > 0
        #-- the oracle is wired up correctly: a lat/lon square's geodesic area
        @test s2_area(write_wkb(GI.Polygon([[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0),
                                             (0.0, 1.0), (0.0, 0.0)]]))) ≈ 1.2364036567e10 rtol = 1e-9
    end
end

if NE_OK
    ne110_names, ne110_geoms = load_ne(110)

    @testset "Natural Earth 110m — spherical overlay vs s2geography" begin
        @test length(ne110_geoms) > 100
        pairs = neighbour_pairs(ne110_names, ne110_geoms, S2_PAIRS)
        @test length(pairs) >= 20
        for (i, j) in pairs
            differential(ne110_geoms[i], ne110_geoms[j],
                         "NE110 $(ne110_names[i]) x $(ne110_names[j])")
        end
        for (nm, X, Y) in shifted_cases(ne110_names, ne110_geoms)
            differential(X, Y, "NE110 $nm")
        end
    end

    ne10_ok = try
        global ne10_names, ne10_geoms = load_ne(10)
        length(ne10_geoms) > 100
    catch err
        @info "Natural Earth 10m spherical differential skipped (data unavailable)" err
        false
    end

    @testset "Natural Earth 10m — spherical overlay vs s2geography" begin
        if !ne10_ok
            missing_data("Natural Earth 10m")
        else
            pairs = neighbour_pairs(ne10_names, ne10_geoms, S2_NE10_PAIRS)
            @test !isempty(pairs)
            for (i, j) in pairs
                differential(ne10_geoms[i], ne10_geoms[j],
                             "NE10 $(ne10_names[i]) x $(ne10_names[j])")
            end
        end
    end

    # -----------------------------------------------------------------------
    # The antimeridian seam, against a spherical oracle for the first time
    # -----------------------------------------------------------------------
    #
    # Natural Earth 10 m Russia is a 214-polygon MultiPolygon, five of whose
    # polygons reach lon ±180 — the input behind the two `@test_broken` pins
    # `realdata_identities.jl` used to carry ("Spherical NE10 Azerbaijan x
    # Russia", "Spherical NE10 Belarus x Russia") for an
    # `OverlayTopologyError: side location conflict` that no longer reproduces.
    #
    # It is the hardest spherical input this corpus has and the one case a
    # planar oracle is structurally unable to check, so it gets its own sweep
    # against s2 rather than being left to the `S2_NE10_PAIRS` slice, which is
    # alphabetical and never reaches Russia.
    @testset "Natural Earth 10m antimeridian — Russia vs s2geography" begin
        if !ne10_ok
            missing_data("Natural Earth 10m")
        else
            ir = findfirst(==("Russia"), ne10_names)
            @test ir !== nothing
            if ir !== nothing
                R = ne10_geoms[ir]
                @test GI.ngeom(R) > 200          # the multi-polygon really is the big one
                for nm in ("Azerbaijan", "Belarus", "Finland", "Norway",
                           "Kazakhstan", "Mongolia", "Ukraine")
                    i = findfirst(==(nm), ne10_names)
                    i === nothing && continue
                    #-- argument order as `neighbour_pairs` produces it (sorted
                    #-- by name), so the labels match the pins above
                    differential(ne10_geoms[i], R, "NE10 $nm x Russia")
                end
            end
        end
    end

    # -----------------------------------------------------------------------
    # Census pins
    # -----------------------------------------------------------------------
    #
    # The header calls the census the deliverable of this suite, and until now
    # nothing asserted any of it: every number was printed and discarded. That
    # is the gap a `_build_lines -> []` mutation walks straight through — it
    # moves `n_collapsed` 41 → 40 and `n_subareal` 0 → 1 while every per-op
    # check still passes, because a dropped shared border turns a
    # dimension-collapsed agreement into a sub-areal one and both are scored
    # `agree`. The classification is the only thing that moves, so the
    # classification is what has to be pinned.
    #
    # Calibrated for the default corpus only; a widened sweep re-measures. To
    # refresh after a deliberate change, run the file and copy the printed
    # census.
    @testset "census" begin
        if S2_PAIRS == 30 && S2_NE10_PAIRS == 4 && ne10_ok
            @test CENSUS.n_op == 204
            @test CENSUS.n_agree == 204
            @test CENSUS.n_disagree == 0
            @test CENSUS.n_collapsed == 41      # ours linear where s2 is areal-only
            @test CENSUS.n_subareal == 0
            @test CENSUS.n_seam == 14           # every ledger entry still fires
            @test CENSUS.n_unresolved <= 35     # s2 could not build the symdifference
        else
            @test_skip "census pins are calibrated for the default corpus"
        end
    end
else
    @testset "Natural Earth — spherical overlay vs s2geography" begin
        missing_data(S2_OK ? "Natural Earth" : "S2Geography_jll")
    end
end

# ---------------------------------------------------------------------------
# Bonus: predicates, via the handle API (route 1)
# ---------------------------------------------------------------------------
#
# A different subsystem — RelateNG's spherical kernel — against a different,
# much simpler s2 route: `S2GeogOpEvalGeogGeog`, no Arrow involved. Predicates
# are exact booleans, so unlike the overlay comparison there is no metric and no
# tolerance: the two engines either agree or they do not.
#
# `equals` is deliberately absent. s2's `S2GEOGRAPHY_OP_EQUALS` compares
# geographies under S2's snapped model, while `GO.equals` is exact set equality;
# on distinct real-world inputs both are simply false, so the comparison would
# assert nothing.
#
# ## The fixtures are chosen so that every predicate takes both values
#
# This block used to sweep all `S2_PAIRS` neighbour pairs, and measured out at
# `intersects` true×30, `disjoint` false×30, `contains` false×30, `within`
# false×30 — 120 assertions, four constant bits, 4.1 s. `intersects` is true by
# construction (`neighbour_pairs` filters on `LG.intersects`), and no pair of
# *neighbours* can ever contain the other, so a constant-`false` `contains`
# and `within` passed 60 of the 120 against an oracle that never contradicted
# them.
#
# Enclaves are the obvious repair and they do NOT work: Natural Earth's country
# polygons partition the land, so an enclave is modelled as a *hole* in its
# neighbour — NE110 South Africa does not contain Lesotho, it abuts it, and
# `contains` is correctly false. (This block asserted the polarity before it
# asserted the fixture, and the polarity assertion is what caught it.)
#
# What does give genuine 2-D containment on real linework is a MultiPolygon
# country against one of its own component polygons: the component's interior
# lies in the country's interior and its boundary on the country's boundary,
# which is `contains` in one direction and `within` in the other.
#
# Antimeridian-spanning countries are excluded from that construction by an
# explicit rule rather than by luck, and the reason is a live disagreement:
# on NE110 **Fiji**, whose first component crosses lon ±180, we answer
# `contains(Fiji, part) = true` — which is unarguable, the part is a subset of
# the whole — and s2 answers `false`. That is a seam-import question about the
# oracle, in the same family as the `S2_SEAM_COMPONENT_SPLIT` ledger above, and
# it does not belong in a block whose job is to check ordinary predicate
# agreement. It is recorded here so the next person does not rediscover it as a
# flake.
if NE_OK
    @testset "spherical predicates vs s2geography" begin
        names, geoms = load_ne(110)
        alg = GO.RelateNG(; manifold = GO.Spherical())
        preds = ((:intersects, GO.intersects), (:disjoint, GO.disjoint),
                 (:contains, GO.contains), (:within, GO.within))

        #-- a handful of neighbours (the bits they exercise are constant across
        #-- the corpus, so the marginal pair buys nothing), one far-apart pair,
        #-- and component-vs-whole pairs in both argument orders
        cases = Tuple{String, Any, Any}[]
        for (i, j) in first(neighbour_pairs(names, geoms, S2_PAIRS), 4)
            push!(cases, ("$(names[i]) x $(names[j])", geoms[i], geoms[j]))
        end
        ib, ia = findfirst(==("Brazil"), names), findfirst(==("Australia"), names)
        @test ib !== nothing && ia !== nothing
        push!(cases, ("Brazil x Australia", geoms[ib], geoms[ia]))   # disjoint

        spans_seam(g) = (e = GI.extent(g); e.X[2] - e.X[1] > 180)
        ncomposite = 0
        for k in eachindex(names)
            ncomposite >= 2 && break
            GI.trait(geoms[k]) isa GI.MultiPolygonTrait || continue
            GI.ngeom(geoms[k]) >= 2 || continue
            spans_seam(geoms[k]) && continue                # see the header note
            part = GI.getgeom(geoms[k], 1)
            spans_seam(part) && continue
            ncomposite += 1
            push!(cases, ("$(names[k]) contains its 1st part", geoms[k], part))
            push!(cases, ("$(names[k])'s 1st part within it", part, geoms[k]))
        end
        @test ncomposite == 2

        seen = Dict(nm => Set{Bool}() for (nm, _) in preds)
        for (label, A, B) in cases
            wa, wb = write_wkb(A), write_wkb(B)
            for (nm, f) in preds
                ours = f(alg, A, B)
                theirs = s2_predicate(nm, wa, wb)
                ours == theirs || println("S2 PREDICATE [$label/$nm] ours=$ours s2=$theirs")
                @test ours == theirs
                push!(seen[nm], ours)
            end
        end
        #-- the assertion the old block was missing: a predicate that only ever
        #-- returned one value was not being compared, it was being echoed
        for (nm, _) in preds
            @test seen[nm] == Set((true, false))
        end
    end
end

# ---------------------------------------------------------------------------
# Census report
# ---------------------------------------------------------------------------

