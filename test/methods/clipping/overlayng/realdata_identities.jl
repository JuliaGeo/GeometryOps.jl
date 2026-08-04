# Oracle-free overlay identities on real data (phase-3 validation, design §4).
#
# No reference engine is involved: every assertion here is an algebraic identity
# the four ops must satisfy among themselves, so the sweep runs unchanged on
# BOTH manifolds — `Planar()` and `Spherical()` — which is the point, since
# there is no spherical oracle to differentiate against.
#
# For each geometry pair (A, B):
#
#   1. conservation   area(A∪B) + area(A∩B) == area(A) + area(B)
#   2. reconstruction A == (A∩B) ∪ (A∖B)
#   3. symdifference  area(A△B) == area(A∖B) + area(B∖A)
#                     area(A△B) == area(A∪B) − area(A∩B)
#   4. validity       every result is `isValid` (planar; the spherical results
#                     get the coordinate-sanity + non-negative-area check that
#                     is meaningful without a spherical validity oracle)
#   5. dimension      on the plane only, each result's per-dimension signature
#                     agrees with GEOS's. Identities 1-3 are all areal, so a
#                     result that silently drops its 1-D components (a shared
#                     border) or its 0-D ones (a corner touch) satisfies every
#                     one of them — neighbouring countries meet along borders,
#                     so that is the common case here, not an exotic one.
#
# Identity 1 is a machine-precision gate since the signed-area fix `69e416484`:
# these are lon/lat coordinates (|x| ≤ 180), so the shoelace has no
# large-coordinate cancellation and `rtol = 1e-12` is a tight bar, not a
# nominal one. (On projected data with coordinates ~1e6 it would NOT be —
# `GO.area(Planar())` runs an untranslated shoelace, unlike JTS/GEOS, and loses
# ~8 digits there; see the note in `xml_suite.jl`.)
#
# Data is availability-gated exactly like `test/wkb.jl`'s Natural Earth block:
# a missing dataset logs and skips, it never fails.
#
# Env knobs (all optional):
#   GO_OVERLAYNG_NE_PAIRS   neighbour pairs per resolution (default 60)
#   GO_OVERLAYNG_NE10       set to 1 to widen the 10 m sweep to the full budget
#   GO_OVERLAYNG_SPH_PAIRS  spherical neighbour pairs (default 8; spherical
#                           overlay is ~5x the cost of planar, and
#                           `s2_differential.jl` sweeps the same 30 pairs
#                           against a real oracle — this file's spherical run
#                           exists to keep *some* spherical real-data coverage
#                           on a platform without `S2Geography_jll`)
#   GO_REQUIRE_DATA         set to 1 to turn a missing dataset into a failure
#                           instead of a skip

using Test
include(joinpath(@__DIR__, "common.jl"))

const NE_PAIRS = something(tryparse(Int, get(ENV, "GO_OVERLAYNG_NE_PAIRS", "")), 60)
const SPH_PAIRS = something(tryparse(Int, get(ENV, "GO_OVERLAYNG_SPH_PAIRS", "")), 8)
const NE10_WIDE = get(ENV, "GO_OVERLAYNG_NE10", "") == "1"

# A missing optional dataset records a skip by default; under `GO_REQUIRE_DATA=1`
# (which CI sets) it fails instead. Without this gate a runner that has lost its
# Natural Earth cache reports the whole file as passing having run zero engine
# assertions — which is how the stale `@test_broken` pins described below
# survived for as long as they did.
const REQUIRE_DATA = get(ENV, "GO_REQUIRE_DATA", "") == "1"
function missing_data(what)
    if REQUIRE_DATA
        println("REQUIRED DATASET MISSING: ", what)
        @test false
    else
        @test_skip what
    end
end

ovl(m, op, a, b) = GO._overlay_ng(m, op, a, b; exact = EX)

function all_finite(g)
    t = GI.trait(g)
    t isa GI.PointTrait && return isfinite(GI.x(g)) && isfinite(GI.y(g))
    t isa GI.GeometryCollectionTrait && return all(all_finite, GI.getgeom(g))
    return all(p -> isfinite(GI.x(p)) && isfinite(GI.y(p)), GI.getpoint(g))
end

# Whether a result can be fed back into the engine as an input on this branch:
# line and area geometries can, points and mixed collections cannot yet.
function reconstructible(g)
    t = GI.trait(g)
    return t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait ||
           t isa GI.LineStringTrait || t isa GI.LinearRingTrait || t isa GI.MultiLineStringTrait
end

# Worst observed conservation residual per manifold, reported at the end.
const WORST = Dict{String, Tuple{Float64, String}}()
function note_worst!(key, rel, label)
    cur = get(WORST, key, (0.0, ""))
    rel > cur[1] && (WORST[key] = (rel, label))
    return nothing
end

"""
    identity_sweep(m, A, B, label; rtol, results)

Run the four identities on one pair. Returns `(ok, detail)`; never throws
(engine errors are reported as a failed check so a single bad pair cannot abort
the sweep).
"""
function identity_sweep(m, A, B, label; rtol = 1e-12)
    try
        I  = ovl(m, GO.OVERLAY_INTERSECTION, A, B)
        U  = ovl(m, GO.OVERLAY_UNION, A, B)
        Dab = ovl(m, GO.OVERLAY_DIFFERENCE, A, B)
        Dba = ovl(m, GO.OVERLAY_DIFFERENCE, B, A)
        SD = ovl(m, GO.OVERLAY_SYMDIFFERENCE, A, B)
        aA, aB = GO.area(m, A), GO.area(m, B)
        aI, aU = GO.area(m, I), GO.area(m, U)
        aDab, aDba, aSD = GO.area(m, Dab), GO.area(m, Dba), GO.area(m, SD)
        scale = aA + aB
        rel(x) = scale == 0 ? abs(x) : abs(x) / scale
        #-- 1. conservation
        r1 = rel(aU + aI - (aA + aB))
        note_worst!(string(typeof(m).name.name), r1, label)
        r1 <= rtol || return (false, "conservation residual $(rel(aU + aI - (aA + aB))) at $label")
        #-- 2. reconstruction: A == (A∩B) ∪ (A∖B), checked by area and, on the
        #--    plane, by GEOS `equals` on the rebuilt geometry — with the same
        #--    graded fallback the differential fuzz uses, because the rebuild
        #--    re-nodes A's edges against emitted (rounded) crossing vertices, so
        #--    an ulp-scale boundary difference from A is expected and benign.
        #--    When A and B meet at a point the intersection is a POINT or a
        #--    mixed collection, which the engine cannot yet take back as an
        #--    *input* (point/GC inputs are phase 3, on a separate branch); the
        #--    identity is then asserted in its area form instead.
        if !reconstructible(I)
            rel(aI + aDab - aA) <= rtol ||
                return (false, "area reconstruction residual $(rel(aI + aDab - aA)) at $label")
            @goto symdiff_checks
        end
        R = ovl(m, GO.OVERLAY_UNION, I, Dab)
        rel(GO.area(m, R) - aA) <= rtol ||
            return (false, "reconstruction area residual $(rel(GO.area(m, R) - aA)) at $label")
        if m isa GO.Planar
            lr, la = result_to_lg(R), result_to_lg(A)
            if !LG.equals(lr, la)
                sd = LG.area(LG.symmetricDifference(lr, la))
                band = 8 * (LG.geomLength(lr) + LG.geomLength(la)) * 180 * eps(Float64)
                note_worst!("reconstruction symdiff / band", sd / band, label)
                sd <= band || return (false,
                    "reconstruction differs from A by symdiff area $sd > rounding band $band at $label")
            end
        end
        #-- 3. symmetric difference
        @label symdiff_checks
        rel(aSD - (aDab + aDba)) <= rtol ||
            return (false, "symdiff-as-two-differences residual $(rel(aSD - (aDab + aDba))) at $label")
        rel(aSD - (aU - aI)) <= rtol ||
            return (false, "symdiff-as-union-minus-intersection residual $(rel(aSD - (aU - aI))) at $label")
        #-- 4. validity, and 5. dimension
        for (nm, g, geos) in (("intersection", I, () -> geos_op(:intersection, A, B)),
                              ("union", U, () -> geos_op(:union, A, B)),
                              ("difference", Dab, () -> geos_op(:difference, A, B)),
                              ("reverse difference", Dba, () -> geos_op(:difference, B, A)),
                              ("symdifference", SD, () -> geos_op(:symdifference, A, B)))
            all_finite(g) || return (false, "non-finite coordinate in $nm at $label")
            if m isa GO.Planar
                lg = result_to_lg(g)
                LG.isValid(lg) ||
                    return (false, "invalid $nm result at $label: " *
                        first(LG.isValidReason(lg), 80))
                #-- the only check in this file that can see a dropped line or
                #-- point: everything above reduces the result to an area, and
                #-- neighbouring countries meet along shared borders, so 1-D
                #-- intersection content is the common case here
                signatures_agree(m, lg, geos()) ||
                    return (false, "$nm dimension signature $(dim_signature(m, lg)) != " *
                        "GEOS $(dim_signature(m, geos())) at $label")
            else
                #-- `GO.area(Planar())` is `abs`-wrapped, so non-negativity only
                #-- has teeth on the sphere, where a ring orientation surviving
                #-- into the result can genuinely come back negative
                GO.area(m, g) >= -1e-12 || return (false, "negative area from $nm at $label")
            end
        end
        return (true, "")
    catch err
        return (false, "ERROR at $label: " * first(sprint(showerror, err), 160))
    end
end

const NE_AVAILABLE = try
    import NaturalEarth, GeoJSON
    true
catch err
    @info "Natural Earth overlay identity sweeps skipped (data unavailable)" err
    false
end

# `load_ne`, `neighbour_pairs`, `shifted` and `shifted_cases` are shared with the
# s2geography differential suite, which sweeps the same corpus.
include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))

# Nothing is pinned: every pair of every sweep must satisfy every identity on
# both manifolds. Two pairs used to be — `Spherical NE10 Azerbaijan x Russia`
# and `Spherical NE10 Belarus x Russia`, for the antimeridian seam defect
# (`OverlayTopologyError: side location conflict` from the 5 of Russia's 214
# polygons that reach lon ±180) — and both cleared long before anyone noticed,
# because they only ever ran in the widened sweep (GO_OVERLAYNG_NE10=1 with a
# large GO_OVERLAYNG_NE_PAIRS) and a `@test_broken` that never executes is never
# reported as unexpectedly passing. That is the argument for `GO_REQUIRE_DATA`
# above: a pin outside the default sweep is a pin nobody is reading.
function run_sweep(m, label, cases; rtol = 1e-12)
    mname = string(typeof(m).name.name)
    for (name, A, B) in cases
        key = "$mname $label $name"          # manifold-specific: a pair can be clean
        ok, detail = identity_sweep(m, A, B, key; rtol)   # on one manifold and not the other
        ok || println("REALDATA FAILURE: ", detail)
        @test ok
    end
    return nothing
end

if NE_AVAILABLE
    ne110_names, ne110_geoms = load_ne(110)

    @testset "Natural Earth 110m — overlay identities" begin
        @test length(ne110_geoms) > 100
        @test !("Sudan" in ne110_names)          # the documented bad input is excluded

        pairs = neighbour_pairs(ne110_names, ne110_geoms, NE_PAIRS)
        @test length(pairs) >= 20
        nbr = [("$(ne110_names[i]) x $(ne110_names[j])", ne110_geoms[i], ne110_geoms[j])
               for (i, j) in pairs]

        #-- shifted-self cases, on countries away from the antimeridian (whose
        #-- overlay is a separate concern handled by `antimeridian_split`)
        shifts = shifted_cases(ne110_names, ne110_geoms)
        @test length(shifts) >= 6

        @testset "planar" begin
            run_sweep(GO.Planar(), "NE110", nbr)
            run_sweep(GO.Planar(), "NE110", shifts)
        end
        @testset "spherical" begin
            run_sweep(GO.Spherical(), "NE110", first(nbr, SPH_PAIRS))
            run_sweep(GO.Spherical(), "NE110", shifts)
        end
    end

    ne10_ok = try
        global ne10_names, ne10_geoms = load_ne(10)
        length(ne10_geoms) > 100
    catch err
        @info "Natural Earth 10m overlay identity sweep skipped (data unavailable)" err
        false
    end

    @testset "Natural Earth 10m — overlay identities" begin
        if !ne10_ok
            missing_data("Natural Earth 10m")
        else
            #-- 10 m rings are one to two orders of magnitude denser, so the CI
            #-- default is a small slice; GO_OVERLAYNG_NE10=1 widens it.
            n = NE10_WIDE ? NE_PAIRS : 8
            pairs = neighbour_pairs(ne10_names, ne10_geoms, n)
            @test !isempty(pairs)
            nbr = [("$(ne10_names[i]) x $(ne10_names[j])", ne10_geoms[i], ne10_geoms[j])
                   for (i, j) in pairs]
            @testset "planar" begin run_sweep(GO.Planar(), "NE10", nbr) end
            @testset "spherical" begin
                run_sweep(GO.Spherical(), "NE10", first(nbr, NE10_WIDE ? SPH_PAIRS : 4))
            end
        end
    end
else
    @testset "Natural Earth — overlay identities" begin missing_data("Natural Earth") end
end

println("Worst overlay conservation residual observed (relative to area(A)+area(B)):")
for (k, (rel, label)) in sort(collect(WORST); by = first)
    println("  $k: $rel   at $label")
end
