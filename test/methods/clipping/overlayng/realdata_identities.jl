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
#   GO_OVERLAYNG_SPH_PAIRS  spherical neighbour pairs (default 30; spherical
#                           overlay is ~5x the cost of planar)

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import Extents

const EX = GO.True()
const NE_PAIRS = something(tryparse(Int, get(ENV, "GO_OVERLAYNG_NE_PAIRS", "")), 60)
const SPH_PAIRS = something(tryparse(Int, get(ENV, "GO_OVERLAYNG_SPH_PAIRS", "")), 30)
const NE10_WIDE = get(ENV, "GO_OVERLAYNG_NE10", "") == "1"

# Natural Earth 110 m **Sudan** is a documented bad input: its exterior ring is
# clockwise, carries a duplicate vertex, and is geodesically self-intersecting.
# s2geography rejects it too. It is excluded by name rather than left to poison
# a sweep — the engine contracts on valid input (design §2.2).
const NE_EXCLUDED = Set(["Sudan"])

ovl(m, op, a, b) = GO._overlay_ng(m, op, a, b; exact = EX)
const OPS = (GO.OVERLAY_INTERSECTION, GO.OVERLAY_UNION, GO.OVERLAY_DIFFERENCE, GO.OVERLAY_SYMDIFFERENCE)

function result_to_lg(g)
    t = GI.trait(g)
    t isa GI.PointTrait && return LG.Point(Float64(GI.x(g)), Float64(GI.y(g)))
    t isa GI.GeometryCollectionTrait &&
        return LG.GeometryCollection(LG.Geometry[result_to_lg(s) for s in GI.getgeom(g)])
    GI.npoint(g) == 0 && return LG.readgeom("GEOMETRYCOLLECTION EMPTY")
    return GI.convert(LG, g)
end

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
        #-- when A and B touch at a single point the intersection is a POINT or a
        #-- mixed collection, which the engine cannot yet take back as an *input*
        #-- (point/GC inputs are phase 3, landing on a separate branch). The
        #-- identity then degenerates to area(A∖B) == area(A), which is asserted
        #-- instead of rebuilding.
        if !reconstructible(I)
            rel(aDab - aA) <= rtol ||
                return (false, "point-touch reconstruction residual $(rel(aDab - aA)) at $label")
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
        #-- 4. validity
        for (nm, g) in (("intersection", I), ("union", U), ("difference", Dab),
                        ("reverse difference", Dba), ("symdifference", SD))
            all_finite(g) || return (false, "non-finite coordinate in $nm at $label")
            GO.area(m, g) >= -1e-12 || return (false, "negative area from $nm at $label")
            if m isa GO.Planar
                LG.isValid(result_to_lg(g)) ||
                    return (false, "invalid $nm result at $label: " *
                        first(LG.isValidReason(result_to_lg(g)), 80))
            end
        end
        return (true, "")
    catch err
        return (false, "ERROR at $label: " * first(sprint(showerror, err), 160))
    end
end

# Pairs whose extents intersect, in a deterministic order (sorted by name), up
# to `limit`. Only pairs that actually share more than a point are kept, so the
# sweep does not fill up with trivially-disjoint neighbours.
function neighbour_pairs(names, geoms, limit)
    order = sortperm(names)
    exts = [GI.extent(g) for g in geoms]
    pairs = Tuple{Int, Int}[]
    for ii in eachindex(order), jj in (ii + 1):length(order)
        i, j = order[ii], order[jj]
        (exts[i] === nothing || exts[j] === nothing) && continue
        Extents.intersects(exts[i], exts[j]) || continue
        push!(pairs, (i, j))
        length(pairs) >= 4 * limit && break
    end
    kept = Tuple{Int, Int}[]
    for (i, j) in pairs
        length(kept) >= limit && break
        try
            LG.intersects(GI.convert(LG, geoms[i]), GI.convert(LG, geoms[j])) || continue
        catch
            continue
        end
        push!(kept, (i, j))
    end
    return kept
end

# Shifted-self pairs: A against a copy translated by a fraction of a degree, so
# every border segment crosses its own copy — the densest crossing workload
# available from this data.
shifted(A, dx, dy) = GO.apply(GI.PointTrait(), A) do p
    (GI.x(p) + dx, GI.y(p) + dy)
end

function load_ne(resolution)
    names = String[]
    geoms = Any[]
    fc = NaturalEarth.naturalearth("admin_0_countries", resolution)
    for f in fc
        g = GeoJSON.geometry(f)
        (g === nothing || GI.npoint(g) == 0) && continue
        nm = try string(f.NAME) catch; "?" end
        nm in NE_EXCLUDED && continue
        t = GI.trait(g)
        (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) || continue
        tg = GO.tuples(g)
        #-- the engine contracts on valid input; NE has a handful of other
        #-- self-intersecting rings besides Sudan, dropped here by the same rule
        try
            LG.isValid(GI.convert(LG, tg)) || continue
        catch
            continue
        end
        push!(names, nm); push!(geoms, tg)
    end
    return names, geoms
end

const NE_AVAILABLE = try
    import NaturalEarth, GeoJSON
    true
catch err
    @info "Natural Earth overlay identity sweeps skipped (data unavailable)" err
    false
end

# Country pairs that the engine currently gets wrong, pinned so that a fix and a
# regression are equally visible. Each is an instance of a defect documented in
# `test/external/jts/overlay_skiplist.jl`.
const NE_KNOWN_DEFECTS = Set{String}([
    #= NE_LEDGER =#
])

function run_sweep(m, label, cases; rtol = 1e-12)
    nfail = 0
    for (name, A, B) in cases
        ok, detail = identity_sweep(m, A, B, "$label $name"; rtol)
        if "$label $name" in NE_KNOWN_DEFECTS
            @test_broken ok
            ok || (nfail += 1)
        else
            ok || println("REALDATA FAILURE: ", detail)
            @test ok
        end
    end
    return nfail
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
        shift_names = ["Brazil", "France", "Egypt", "Australia", "Chile", "Norway",
                       "Indonesia", "India", "Kazakhstan", "Argentina"]
        shifts = Any[]
        for nm in shift_names
            idx = findfirst(==(nm), ne110_names)
            idx === nothing && continue
            A = ne110_geoms[idx]
            push!(shifts, ("$nm shifted", A, shifted(A, 0.5, 0.25)))
        end
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
            @test_skip "Natural Earth 10m unavailable"
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
end

# GADM (full-resolution national boundaries) is only reachable from the docs
# environment — it is not a test dependency, and its GeoPackages download on
# first use. Gated exactly like the Natural Earth block: available or skipped,
# never failing.
@testset "GADM — overlay identities" begin
    gadm_ok = try
        import GADM
        true
    catch err
        @info "GADM overlay identity sweep skipped (GADM.jl unavailable in the test environment)" err
        false
    end
    if !gadm_ok
        @test_skip "GADM unavailable"
    else
        cases = Any[]
        for (x, y) in (("EGY", "SDN"), ("FRA", "ITA"), ("GRC", "TUR"))
            try
                A = GO.tuples(GI.getgeom(GADM.get(x), 1))
                B = GO.tuples(GI.getgeom(GADM.get(y), 1))
                push!(cases, ("$x x $y", A, B))
            catch err
                @info "GADM pair unavailable" x y err
            end
        end
        @test !isempty(cases)
        @testset "planar" begin run_sweep(GO.Planar(), "GADM", cases) end
        @testset "spherical" begin run_sweep(GO.Spherical(), "GADM", cases) end
    end
end

println("Worst overlay conservation residual observed (relative to area(A)+area(B)):")
for (k, (rel, label)) in sort(collect(WORST); by = first)
    println("  $k: $rel   at $label")
end
