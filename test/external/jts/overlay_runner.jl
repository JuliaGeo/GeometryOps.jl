# Runner for the vendored JTS overlay XML conformance suites.
#
# Parameterized over the overlay implementation (like `relate_runner.jl` is over
# relate), so it stays engine-agnostic: it never calls GeometryOps for its own
# pass/fail logic, only through the `overlay_fn` it is handed.
#
# Two entry points:
#
# - `run_overlay_cases`   — expected-result conformance: every `intersection` /
#   `union` / `difference` / `symdifference` op (and their `*NG` spellings) is
#   run and compared against the XML's expected geometry with GEOS topological
#   `equals`. Coordinate equality is deliberately NOT used: this engine emits an
#   exact arrangement rounded once at emission, so its vertex coordinates need
#   not be bit-identical to JTS's.
# - `run_overlay_identity_cases` — oracle-free identity sweep (JTS's own
#   `overlayAreaTest` / `TestCaseGeometryFunctions.areaDelta`): for every case
#   geometry pair, all five overlay area identities must hold. Used on the
#   `robust/overlay` corpus, whose XML files mostly carry no expected results
#   (they are regression *inputs* harvested from GEOS/PostGIS/QGIS/shapely bug
#   reports).

isdefined(@__MODULE__, :load_test_run) || include(joinpath(@__DIR__, "jts_testset_reader.jl"))
isdefined(@__MODULE__, :OVERLAY_SKIPLIST) || include(joinpath(@__DIR__, "overlay_skiplist.jl"))

using Test
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

# XML op name (lowercased) → canonical op symbol. The `*NG` spellings are the
# same operations routed through JTS's OverlayNG rather than the legacy engine;
# both are expected to give the same answer on FLOATING precision.
const OVERLAY_OP_NAMES = Dict{String, Symbol}(
    "intersection"    => :intersection, "intersectionng"  => :intersection,
    "union"           => :union,        "unionng"         => :union,
    "difference"      => :difference,   "differenceng"    => :difference,
    "symdifference"   => :symdifference,"symdifferenceng" => :symdifference,
)

const OVERLAY_OP_ORDER = (:intersection, :union, :difference, :symdifference)

# -- geometry plumbing --------------------------------------------------------

# `jts_wkt_to_geom` returns `GO.tuples` output for ordinary WKT and a raw
# LibGEOS geometry for WKT it cannot handle (EMPTY, GEOMETRYCOLLECTION,
# LINEARRING). LinearRings are converted here so the engine sees a native
# geometry; EMPTY and GEOMETRYCOLLECTION are left alone (empty LibGEOS geoms
# answer `GI.isempty` correctly, which is what the engine's short circuits
# need; collections are skipped).
function jts_overlay_operand(g)
    g === nothing && return nothing
    GI.trait(g) isa GI.LinearRingTrait && return GO.tuples(g)
    return g
end

# Classification used for skip decisions, so an unsupported input is recorded
# with a precise reason instead of surfacing as an engine error.
function overlay_operand_kind(g)
    t = GI.trait(g)
    (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) && return :area
    (t isa GI.LineStringTrait || t isa GI.LinearRingTrait ||
        t isa GI.MultiLineStringTrait) && return :line
    (t isa GI.PointTrait || t isa GI.MultiPointTrait) && return :point
    t isa GI.GeometryCollectionTrait && return :collection
    return :unsupported
end

_is_lg(g) = g isa LG.AbstractGeometry

# Emptiness that works for every representation the corpus and the engine
# produce: LibGEOS geometries, bare coordinate tuples (what `GO.tuples` yields
# for POINT WKT), GeoInterface wrappers, and geometry collections (for which
# `GI.npoint` has no method).
function overlay_is_empty(g)
    _is_lg(g) && return LG.isEmpty(g)
    t = GI.trait(g)
    t isa GI.PointTrait && return false
    t isa GI.GeometryCollectionTrait && return all(overlay_is_empty, GI.getgeom(g))
    return GI.npoint(g) == 0
end

# Geometry → LibGEOS, for the comparison oracle only. Empty geometries are
# mapped to one canonical empty: GEOS `equals` is `true` between any two empty
# geometries regardless of declared type, so empty-result *type* is checked
# separately (`overlay_result_dimension`), not through `equals`.
function overlay_to_lg(g)
    _is_lg(g) && return g
    t = GI.trait(g)
    t isa GI.PointTrait && return LG.Point(Float64(GI.x(g)), Float64(GI.y(g)))
    overlay_is_empty(g) && return LG.readgeom("GEOMETRYCOLLECTION EMPTY")
    t isa GI.GeometryCollectionTrait &&
        return LG.GeometryCollection(LG.Geometry[overlay_to_lg(sub) for sub in GI.getgeom(g)])
    return GI.convert(LG, g)
end

overlay_wkt(g) = try
    LG.writegeom(overlay_to_lg(g))
catch err
    "<unrepresentable: " * first(sprint(showerror, err), 80) * ">"
end

# Topological dimension of a (possibly empty, possibly heterogeneous) geometry:
# the max dimension over its components, or `-1` for an empty geometry with no
# components. Used to check the *type* of empty results, which `equals` cannot.
function overlay_result_dimension(g)
    t = GI.trait(g)
    if t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait
        return 2
    elseif t isa GI.LineStringTrait || t isa GI.LinearRingTrait || t isa GI.MultiLineStringTrait
        return 1
    elseif t isa GI.PointTrait || t isa GI.MultiPointTrait
        return 0
    elseif t isa GI.GeometryCollectionTrait
        d = -1
        for sub in GI.getgeom(g)
            d = max(d, overlay_result_dimension(sub))
        end
        return d
    end
    return -1
end

# GEOS `equals`, with the empty/empty case answered directly (some GEOS
# versions are fussy about `equals` on mixed-type empties).
function overlay_geom_equals(ours, expected)
    eo, ee = overlay_is_empty(ours), overlay_is_empty(expected)
    (eo || ee) && return eo && ee
    return LG.equals(overlay_to_lg(ours), overlay_to_lg(expected))
end

# Diagnostic for a failing comparison: the area of the symmetric difference of
# the two results (0 for a pure boundary/orientation difference, small for a
# rounding-scale difference, large for a genuine divergence).
function overlay_symdiff_area(ours, expected)
    try
        (overlay_is_empty(ours) || overlay_is_empty(expected)) && return NaN
        return LG.area(LG.symmetricDifference(overlay_to_lg(ours), overlay_to_lg(expected)))
    catch
        return NaN
    end
end

# -- expected-result conformance ---------------------------------------------

"""
    run_overlay_cases(overlay_fn, files; skiplist = OVERLAY_SKIPLIST,
                      point_inputs = false, record_tests = true)

Run the overlay ops of the JTS XML files in `files` against an overlay
implementation `overlay_fn(op::Symbol, a, b)` (`op` ∈ `:intersection`,
`:union`, `:difference`, `:symdifference`) and compare each result to the XML's
expected geometry with GEOS topological `equals`.

Ops are skipped, never silently dropped, when:

- their `(file, case_index, op_name, arg_order)` key is in `skiplist`
  (every entry carries a written justification — see `overlay_skiplist.jl`);
- the run declares a non-FLOATING precision model (fixed-precision overlay is
  permanently out of scope: the engine has no precision model by design);
- the op is not one of the four overlay ops;
- an operand is a `GEOMETRYCOLLECTION` (unsupported input);
- an operand is a point geometry and `point_inputs = false`.

`arg_order` is `"AB"` when the op's `arg1` is the case's A geometry and `"BA"`
when it is B (JTS runs `difference` both ways round in the `*NG` files).

Each executed op contributes one `@test` unless `record_tests = false` (used by
the exploratory driver). Returns `(; per_file, skipped, failures)`.
"""
function run_overlay_cases(overlay_fn, files;
        skiplist::Set{Tuple{String, Int, String, String}} = OVERLAY_SKIPLIST,
        point_inputs::Bool = false, record_tests::Bool = true)
    per_file = NamedTuple{(:file, :n_pass, :n_fail, :n_skip), NTuple{4, Any}}[]
    skipped = NamedTuple{(:file, :case_index, :description, :op, :arg_order, :reason),
        NTuple{6, Any}}[]
    failures = NamedTuple{(:file, :case_index, :description, :op, :arg_order,
        :detail, :symdiff_area), NTuple{7, Any}}[]
    for filepath in files
        file = basename(filepath)
        run = load_test_run(filepath)
        n_pass = n_fail = n_skip = 0
        for (case_index, case) in enumerate(run.cases)
            a = jts_overlay_operand(case.geom_a)
            b = jts_overlay_operand(case.geom_b)
            for item in case.items
                op_name = item.operation
                arg_order = item.arg1 === case.geom_a ? "AB" : "BA"
                skip!(reason) = begin
                    n_skip += 1
                    push!(skipped, (; file, case_index, description = case.description,
                        op = op_name, arg_order, reason))
                end
                if (file, case_index, op_name, arg_order) in skiplist
                    skip!("in skiplist"); continue
                elseif run.precision_model != "FLOATING"
                    skip!("non-FLOATING precision model ($(run.precision_model))"); continue
                end
                op = get(OVERLAY_OP_NAMES, lowercase(op_name), nothing)
                if op === nothing
                    skip!("non-overlay op"); continue
                elseif item.arg2 === nothing || b === nothing
                    skip!("binary op without a B geometry"); continue
                end
                x, y = arg_order == "AB" ? (a, b) : (b, a)
                kx, ky = overlay_operand_kind(x), overlay_operand_kind(y)
                if kx === :collection || ky === :collection
                    skip!("geometry collection input unsupported"); continue
                elseif kx === :unsupported || ky === :unsupported
                    skip!("unsupported input geometry"); continue
                elseif (kx === :point || ky === :point) && !point_inputs
                    skip!("point input unsupported"); continue
                end
                passed = false
                detail = ""
                sdarea = NaN
                try
                    ours = overlay_fn(op, x, y)
                    passed = overlay_geom_equals(ours, item.expected_result)
                    if passed
                        #-- `equals` cannot see the type of an empty result, so
                        #-- check its dimension separately (that is exactly what
                        #-- TestOverlayEmpty / TestNGOverlayEmpty are about).
                        if overlay_is_empty(item.expected_result)
                            de = overlay_result_dimension(item.expected_result)
                            do_ = overlay_result_dimension(ours)
                            #-- JTS returns GEOMETRYCOLLECTION EMPTY (dim -1) where
                            #-- the result dimension is not determined; accept that.
                            passed = de < 0 || do_ < 0 || de == do_
                            passed || (detail = "empty result dimension $do_ != expected $de")
                        end
                    else
                        sdarea = overlay_symdiff_area(ours, item.expected_result)
                        detail = string("not equals (symdiff area ", sdarea, ")\n      ours = ",
                            first(overlay_wkt(ours), 400), "\n      want = ",
                            first(overlay_wkt(item.expected_result), 400))
                    end
                catch err
                    passed = false
                    detail = "ERROR " * first(sprint(showerror, err), 200)
                end
                if passed
                    n_pass += 1
                else
                    n_fail += 1
                    push!(failures, (; file, case_index, description = case.description,
                        op = op_name, arg_order, detail, symdiff_area = sdarea))
                end
                record_tests && @test passed
            end
        end
        push!(per_file, (; file, n_pass, n_fail, n_skip))
    end
    return (; per_file, skipped, failures)
end

# -- oracle-free identity sweep (JTS `overlayAreaTest`) -----------------------

"""
    overlay_area_delta(area_fn, overlay_fn, a, b)

Port of JTS `TestCaseGeometryFunctions.areaDelta` (the engine behind the
`overlayAreaTest` op used throughout `robust/overlay`): the maximum absolute
violation of the five overlay area identities

    A       = (A∩B) + (A∖B)
    B       = (A∩B) + (B∖A)
    (A△B)   = (A∖B) + (B∖A)
    (A∪B)   = (A∩B) + (A△B)
    (A∪B)   = (A∩B) + (A∖B) + (B∖A)

Returns `(delta, total)` where `total = area(A) + area(B)` (the natural scale
for the residual). `delta` is `0.0` when either input has zero area, exactly as
JTS does — the identities are vacuous for non-areal inputs.
"""
function overlay_area_delta(area_fn, overlay_fn, a, b)
    areaA = area_fn(a)
    areaB = area_fn(b)
    (areaA == 0 || areaB == 0) && return (0.0, areaA + areaB)
    areaU   = area_fn(overlay_fn(:union, a, b))
    areaI   = area_fn(overlay_fn(:intersection, a, b))
    areaDab = area_fn(overlay_fn(:difference, a, b))
    areaDba = area_fn(overlay_fn(:difference, b, a))
    areaSD  = area_fn(overlay_fn(:symdifference, a, b))
    delta = max(
        abs(areaA - areaI - areaDab),
        abs(areaB - areaI - areaDba),
        abs(areaDab + areaDba - areaSD),
        abs(areaI + areaSD - areaU),
        abs(areaU - areaI - areaDab - areaDba))
    return (delta, areaA + areaB)
end

"""
    run_overlay_identity_cases(overlay_fn, area_fn, files; kwargs...)

Run the oracle-free overlay identity sweep over every `<case>` of `files`,
regardless of which ops the XML declares (the `robust/overlay` corpus is a
collection of hard *inputs*; most of its cases carry no expected geometry).

For each case with two area operands the five identities of
[`overlay_area_delta`](@ref) must hold to `rtol` relative to `area(A)+area(B)`,
and every one of the five overlay results must satisfy `result_valid_fn`. Cases
are skipped, with a recorded reason, when an operand is unsupported, when an
input fails `valid_fn` (the engine contracts on valid input, design §2.2), or
when the case is in `skiplist` (keyed `(file, case_index, "identity", "AB")`).

Returns `(; per_file, skipped, failures, worst)` where `worst` is the largest
observed `(relative delta, file, case_index)`.
"""
function run_overlay_identity_cases(overlay_fn, area_fn, files;
        valid_fn = _ -> true, result_valid_fn = _ -> true,
        skiplist::Set{Tuple{String, Int, String, String}} = OVERLAY_SKIPLIST,
        rtol::Float64 = 1e-10, record_tests::Bool = true, max_points::Int = typemax(Int))
    per_file = NamedTuple{(:file, :n_pass, :n_fail, :n_skip), NTuple{4, Any}}[]
    skipped = NamedTuple{(:file, :case_index, :description, :reason), NTuple{4, Any}}[]
    failures = NamedTuple{(:file, :case_index, :description, :detail), NTuple{4, Any}}[]
    worst = (rel = 0.0, file = "", case_index = 0)
    for filepath in files
        file = basename(filepath)
        run = load_test_run(filepath)
        n_pass = n_fail = n_skip = 0
        for (case_index, case) in enumerate(run.cases)
            skip!(reason) = begin
                n_skip += 1
                push!(skipped, (; file, case_index, description = case.description, reason))
            end
            if (file, case_index, "identity", "AB") in skiplist
                skip!("in skiplist"); continue
            elseif run.precision_model != "FLOATING"
                skip!("non-FLOATING precision model ($(run.precision_model))"); continue
            end
            a = jts_overlay_operand(case.geom_a)
            b = jts_overlay_operand(case.geom_b)
            if b === nothing
                skip!("case has no B geometry"); continue
            end
            ka, kb = overlay_operand_kind(a), overlay_operand_kind(b)
            if ka !== :area || kb !== :area
                skip!("identity sweep needs two area inputs (got $ka/$kb)"); continue
            elseif GI.npoint(a) == 0 || GI.npoint(b) == 0
                skip!("empty input"); continue
            elseif GI.npoint(a) + GI.npoint(b) > max_points
                skip!("over the $max_points-point budget"); continue
            elseif !(valid_fn(a) && valid_fn(b))
                skip!("invalid input (engine contracts on valid input)"); continue
            end
            passed = true
            detail = ""
            try
                invalid = String[]
                checked(o, x, y) = begin
                    r = overlay_fn(o, x, y)
                    result_valid_fn(r) || push!(invalid, string(o))
                    r
                end
                delta, total = overlay_area_delta(area_fn, checked, a, b)
                rel = total == 0 ? 0.0 : delta / total
                rel > worst.rel && (worst = (rel = rel, file = file, case_index = case_index))
                if !(rel <= rtol)
                    passed = false
                    detail = "area identity residual $delta (relative $rel) over total $total"
                elseif !isempty(invalid)
                    passed = false
                    detail = "invalid result geometry from " * join(unique(invalid), ", ")
                end
            catch err
                passed = false
                detail = "ERROR " * first(sprint(showerror, err), 300)
            end
            if passed
                n_pass += 1
            else
                n_fail += 1
                push!(failures, (; file, case_index, description = case.description, detail))
            end
            record_tests && @test passed
        end
        push!(per_file, (; file, n_pass, n_fail, n_skip))
    end
    return (; per_file, skipped, failures, worst)
end
