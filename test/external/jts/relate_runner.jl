# Runner for the vendored JTS relate XML conformance suite.
# Parameterized over the relate implementation so it can be smoke-tested
# before the engine exists (Task 14); fully activated in Task 23.

isdefined(@__MODULE__, :load_test_run) || include(joinpath(@__DIR__, "jts_testset_reader.jl"))
isdefined(@__MODULE__, :RELATE_SKIPLIST) || include(joinpath(@__DIR__, "relate_skiplist.jl"))

using Test

# Standard DE-9IM pattern matching over the 9-character string forms, kept
# local so the runner stays engine-agnostic (it must not depend on the
# implementation under test for its own pass/fail logic).
_de9im_entry_matches(c::Char, p::Char) =
    p == '*' ? true :
    p == 'T' ? c in ('0', '1', '2') :
    c == p
_de9im_matches(im::AbstractString, pattern::AbstractString) =
    all(_de9im_entry_matches(im[i], pattern[i]) for i in 1:9)

"""
    run_relate_cases(relate_fn, pattern_fn, predicate_fns, files; skiplist = RELATE_SKIPLIST, check_matrix = true)

Run the JTS relate XML test cases in `files` (paths to vendored XML files)
against a relate implementation:

- `relate_fn(a, b)::DE9IM` — computes the full DE-9IM matrix. For every
  `relate` op (unless `check_matrix = false`), the full matrix's `string`
  form is matched against the op's pattern as a consistency check: the
  matrix route and the pattern-predicate route must agree on the expected
  result for the op to pass.
- `pattern_fn(a, b, pattern)::Bool` — evaluates a DE-9IM pattern match
  (used for `relate` ops, whose `arg3` is the pattern).
- `predicate_fns::AbstractDict` — maps lowercase op names (`"intersects"`,
  `"contains"`, ...) to `(a, b) -> Bool` closures. Ops with no entry are
  recorded as skipped.

Ops whose `(file, case_index, op, arg_order)` key is in `skiplist` are recorded
as skipped, never silently dropped — as are ops outside `BOOLEAN_OPS` (e.g.
`getboundary`) and runs with a FIXED precision model. `case_index` is the
1-based index of the `<case>` within the file; `arg_order` is `"AB"` when the
op's `arg1` is the case's A geometry and `"BA"` when it is B (some files, e.g.
TestRelateGC.xml, run each predicate both ways, and case descriptions are not
unique within a file, so neither alone can serve as a key).

Each executed op contributes one `@test`. Returns a NamedTuple
`(per_file, skipped)` where `per_file` is a `Vector` of
`(file, n_pass, n_fail, n_skip)` NamedTuples and `skipped` is a `Vector` of
`(file, case_index, description, op, arg_order, reason)` NamedTuples for every
skipped op.
"""
function run_relate_cases(relate_fn, pattern_fn, predicate_fns, files;
        skiplist::Set{Tuple{String, Int, String, String}} = RELATE_SKIPLIST,
        check_matrix::Bool = true)
    per_file = NamedTuple{(:file, :n_pass, :n_fail, :n_skip), Tuple{String, Int, Int, Int}}[]
    skipped = NamedTuple{(:file, :case_index, :description, :op, :arg_order, :reason),
        Tuple{String, Int, String, String, String, String}}[]
    for filepath in files
        file = basename(filepath)
        run = load_test_run(filepath)
        n_pass = 0
        n_fail = 0
        n_skip = 0
        skip!(case_index, case, item, arg_order, reason) = begin
            n_skip += 1
            push!(skipped, (; file, case_index, description = case.description,
                op = item.operation, arg_order, reason))
        end
        @testset "$file" begin
            for (case_index, case) in enumerate(run.cases)
                @testset "$(case.description)" begin
                    for item in case.items
                        op = lowercase(item.operation)
                        # "AB" when the op's arg1 is the case's A geometry, "BA"
                        # when it is B (e.g. TestRelateGC runs each predicate
                        # with the arguments both ways round).
                        arg_order = item.arg1 === case.geom_a ? "AB" : "BA"
                        if (file, case_index, item.operation, arg_order) in skiplist
                            skip!(case_index, case, item, arg_order, "in skiplist")
                            continue
                        elseif run.precision_model != "FLOATING"
                            skip!(case_index, case, item, arg_order, "non-FLOATING precision model ($(run.precision_model))")
                            continue
                        elseif !(op in BOOLEAN_OPS)
                            skip!(case_index, case, item, arg_order, "non-boolean op")
                            continue
                        end
                        a, b = item.arg1, item.arg2
                        passed = false
                        try
                            if op == "relate"
                                passed = pattern_fn(a, b, item.pattern) == item.expected_result
                                if passed && check_matrix
                                    # Consistency check: the full matrix, matched
                                    # against the pattern, must reproduce the
                                    # pattern-predicate result.
                                    im_str = string(relate_fn(a, b))
                                    passed = _de9im_matches(im_str, item.pattern) == item.expected_result
                                end
                            elseif haskey(predicate_fns, op)
                                passed = predicate_fns[op](a, b) == item.expected_result
                            else
                                skip!(case_index, case, item, arg_order, "no predicate function provided")
                                continue
                            end
                        catch err
                            @error "relate case errored" file case.description item.operation exception = (err, catch_backtrace())
                        end
                        passed ? (n_pass += 1) : (n_fail += 1)
                        @test passed
                    end
                end
            end
        end
        push!(per_file, (; file, n_pass, n_fail, n_skip))
    end
    return (; per_file, skipped)
end
