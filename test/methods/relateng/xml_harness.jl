using Test

include(joinpath(@__DIR__, "..", "..", "external", "jts", "jts_testset_reader.jl"))
include(joinpath(@__DIR__, "..", "..", "external", "jts", "relate_runner.jl"))

const JTS_DATA_DIR = joinpath(@__DIR__, "..", "..", "data", "jts")

@testset "relate XML parsing" begin
    cases = load_test_cases(joinpath(JTS_DATA_DIR, "general", "TestRelatePP.xml"))
    @test !isempty(cases)
    item = first(first(cases).items)
    @test item.operation == "relate"
    @test item.expected_result isa Bool       # boolean ops parse as Bool now
    @test item.pattern isa String && length(item.pattern) == 9
end

@testset "run-level metadata" begin
    # validate files declare <precisionModel type="FLOATING"/>
    run = load_test_run(joinpath(JTS_DATA_DIR, "validate", "TestRelateAA.xml"))
    @test run.precision_model == "FLOATING"
    @test run.n_skipped_ops == 0
    @test !isempty(run.cases)

    # robust files declare a scale => FIXED precision model
    run_fixed = load_test_run(joinpath(JTS_DATA_DIR, "validate", "TestRobustRelate.xml"))
    @test run_fixed.precision_model == "FIXED"

    # misc files have a run-level <desc> and no precisionModel (=> FLOATING default)
    run_desc = load_test_run(joinpath(JTS_DATA_DIR, "general", "TestRelateEmpty.xml"))
    @test run_desc.precision_model == "FLOATING"
    @test run_desc.description isa String && !isempty(run_desc.description)
end

@testset "cases without <b> and non-boolean ops" begin
    # TestBoundary cases have only an <a> geometry and a unary getboundary op
    # whose expected result is a geometry.
    cases = load_test_cases(joinpath(JTS_DATA_DIR, "general", "TestBoundary.xml"))
    @test !isempty(cases)
    case = first(cases)
    @test case.geom_b === nothing
    item = first(case.items)
    @test item.operation == "getboundary"
    @test item.arg2 === nothing
    @test item.pattern === nothing
    @test !(item.expected_result isa Bool)   # parses as a geometry
end

@testset "every vendored relate file parses" begin
    for dir in ("general", "validate"), file in readdir(joinpath(JTS_DATA_DIR, dir))
        endswith(file, ".xml") || continue
        run = load_test_run(joinpath(JTS_DATA_DIR, dir, file))
        @test !isempty(run.cases)
        @test run.n_skipped_ops == 0
        @test all(c -> !isempty(c.items), run.cases)
    end
end

@testset "relate runner skiplist machinery" begin
    # Smoke-test the runner shape with every op skipped and always-throwing
    # closures: nothing should execute, everything should be recorded.
    file = joinpath(JTS_DATA_DIR, "general", "TestRelatePP.xml")
    cases = load_test_cases(file)
    skiplist = Set(
        (basename(file), case_index, i.operation, i.arg1 === c.geom_a ? "AB" : "BA")
            for (case_index, c) in enumerate(cases) for i in c.items)
    boom(args...) = error("must not be called while everything is skipped")
    summary = run_relate_cases(boom, boom, Dict{String, Function}(), [file]; skiplist)
    @test length(summary.per_file) == 1
    stats = only(summary.per_file)
    @test stats.file == basename(file)
    @test stats.n_pass == 0 && stats.n_fail == 0
    @test stats.n_skip == sum(c -> length(c.items), cases)
    @test length(summary.skipped) == stats.n_skip
    @test all(s -> s.reason == "in skiplist", summary.skipped)
end

@testset "relate runner execute path" begin
    # Exercise the runner's execute path with stub implementations (not a real
    # engine, so we only assert count consistency, not specific pass counts —
    # except that TestRelatePP's ops all expect `true`, which the stubs return,
    # so no @test inside the runner fails and pollutes this testset).
    file = joinpath(JTS_DATA_DIR, "general", "TestRelatePP.xml")
    cases = load_test_cases(file)
    n_ops = sum(c -> length(c.items), cases)
    relate_fn = (a, b) -> error("relate_fn is unused until Task 23")
    pattern_fn = (a, b, p) -> true
    predicate_fns = Dict{String, Function}(op => ((a, b) -> true) for op in BOOLEAN_OPS if op != "relate")
    summary = run_relate_cases(relate_fn, pattern_fn, predicate_fns, [file])
    stats = only(summary.per_file)
    @test stats.file == basename(file)
    # nothing skipped: FLOATING precision, all ops boolean, all predicates provided
    @test stats.n_skip == 0
    @test isempty(summary.skipped)
    # every boolean op executed exactly once...
    @test stats.n_pass + stats.n_fail == n_ops
    # ...so every op of every case was visited (executed or skipped)
    @test stats.n_pass + stats.n_fail + stats.n_skip == n_ops
end
