# Full JTS relate XML conformance suite (Task 23): every vendored case in
# test/data/jts/{general,validate} run against the real RelateNG engine.
# The parser/runner machinery itself is unit-tested in xml_harness.jl; this
# file is the conformance run.
#
# Every executed op contributes one `@test` inside `run_relate_cases`. On top
# of that, the expected pass/skip counts per file are pinned below so that an
# accidental mass-skip (e.g. a parser regression silently dropping ops) fails
# loudly instead of shrinking the suite.
#
# Runtime: the suite itself is fast (~0.5 s warm for all 6458 ops); the
# dominant cost is one-time compilation of the engine for the LibGEOS-backed
# geometry types used by EMPTY/GC cases (~10 min cold on an M-series laptop,
# paid once per test process and shared with the other relateng test files).

using Test
import GeometryOps as GO

include(joinpath(@__DIR__, "..", "..", "external", "jts", "jts_testset_reader.jl"))
include(joinpath(@__DIR__, "..", "..", "external", "jts", "relate_runner.jl"))

const JTS_DATA_DIR = joinpath(@__DIR__, "..", "..", "data", "jts")

const RELATENG = GO.RelateNG()

relateng_relate(a, b) = GO.relate(RELATENG, a, b)
relateng_pattern(a, b, pattern) = GO.relate(RELATENG, a, b, pattern)

# Predicate factories return fresh mutable predicate state, so each call
# constructs its own predicate. `equals` (used by some JTS files) is JTS
# `equalsTopo`.
const RELATENG_PREDICATE_FACTORIES = Dict{String, Function}(
    "intersects" => GO.pred_intersects,
    "disjoint"   => GO.pred_disjoint,
    "contains"   => GO.pred_contains,
    "within"     => GO.pred_within,
    "covers"     => GO.pred_covers,
    "coveredby"  => GO.pred_coveredby,
    "crosses"    => GO.pred_crosses,
    "touches"    => GO.pred_touches,
    "overlaps"   => GO.pred_overlaps,
    "equalstopo" => GO.pred_equalstopo,
    "equals"     => GO.pred_equalstopo,
)

relateng_predicate_fns() = Dict{String, Function}(
    name => ((a, b) -> GO.relate_predicate(RELATENG, factory(), a, b))
    for (name, factory) in RELATENG_PREDICATE_FACTORIES)

_xml_files(dir) = sort!(filter!(f -> endswith(f, ".xml"),
    readdir(joinpath(JTS_DATA_DIR, dir); join = true)))

# Expected (pass, skip) counts per file, pinned as of Task 23.
# Skips are: TestBoundary's 12 unary geometry-valued `getboundary` ops
# (not relate ops at all — the file has no <op name="relate"> with a boundary
# rule, just boundary construction cases); TestRobustRelate's single op (the
# run declares a FIXED precision model, which assumes snapped coordinates);
# and the 2 documented skiplist entries (see test/external/jts/relate_skiplist.jl).
const EXPECTED_COUNTS = [
    ("general", "TestBoundary.xml",           0, 12),
    ("general", "TestRelateAA.xml",          41,  0),
    ("general", "TestRelateEmpty.xml",      572,  0),
    ("general", "TestRelateGC.xml",         328,  0),
    ("general", "TestRelateLA.xml",          13,  0),
    ("general", "TestRelateLL.xml",          45,  0),
    ("general", "TestRelatePA.xml",         121,  0),
    ("general", "TestRelatePL.xml",           8,  0),
    ("general", "TestRelatePP.xml",           4,  0),
    ("validate", "TestRelateAA-big.xml",      2,  0),
    ("validate", "TestRelateAA.xml",       1177,  0),
    ("validate", "TestRelateAC.xml",         11,  0),
    ("validate", "TestRelateLA.xml",        847,  0),
    ("validate", "TestRelateLC.xml",         22,  0),
    ("validate", "TestRelateLL.xml",       1584,  0),
    ("validate", "TestRelatePA.xml",        451,  0),
    ("validate", "TestRelatePL.xml",       1088,  1),
    ("validate", "TestRelatePP.xml",        143,  0),
    ("validate", "TestRobustRelate.xml",      0,  1),
    ("validate", "TestRobustRelateFloat.xml", 1,  1),
]

@testset "JTS relate XML conformance suite" begin
    for dir in ("general", "validate")
        files = _xml_files(dir)
        expected = [(f, np, ns) for (d, f, np, ns) in EXPECTED_COUNTS if d == dir]
        @test basename.(files) == first.(expected)
        @testset "$dir" begin
            summary = run_relate_cases(relateng_relate, relateng_pattern,
                relateng_predicate_fns(), files)
            @test sum(s -> s.n_fail, summary.per_file) == 0
            for (s, (file, n_pass, n_skip)) in zip(summary.per_file, expected)
                @test s.file == file
                @test (s.file, s.n_pass) == (file, n_pass)
                @test (s.file, s.n_skip) == (file, n_skip)
            end
            # Every skip is one of the three accounted-for kinds.
            for sk in summary.skipped
                @test sk.reason in ("in skiplist", "non-boolean op",
                    "non-FLOATING precision model (FIXED)")
            end
        end
    end
end
