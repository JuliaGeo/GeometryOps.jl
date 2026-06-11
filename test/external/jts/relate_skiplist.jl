# Skiplist for the JTS relate XML conformance suite (`relate_runner.jl`).
#
# Each entry is a `(file, case_index, op name, arg order)` tuple identifying one
# `<op>` in one `<case>` of one vendored XML file:
#
# - `file` is the basename, e.g. "TestRelateAA.xml".
# - `case_index` is the 1-based index of the `<case>` within the file.
# - `op name` is as written in the XML, e.g. "equalsTopo".
# - `arg order` is "AB" when the op's `arg1` is the case's A geometry and "BA"
#   when it is B (e.g. TestRelateGC.xml runs each predicate with the arguments
#   both ways round).
#
# Case descriptions are NOT part of the key — they are not unique within a file
# (e.g. TestRelateLL.xml has two distinct cases both described as "Line vs line
# - pointwise equal") — so the human-readable description belongs in the
# justification comment instead.
#
# DISCIPLINE: every entry MUST be accompanied by a comment giving the case
# description and explaining exactly why GeometryOps diverges from JTS on that
# case (or why the case cannot run), ideally with an issue/task reference.
# Entries without a justification comment must not be merged. Skipped cases are
# reported by `run_relate_cases` — they are never silently dropped.
#
# Example entry:
#     # "P/L - empty point VS empty line" (case 3): GO returns FF* for
#     # empty-geometry boundary, JTS expects F0*; tracked in #XXXX
#     ("TestRelateEmpty.xml", 3, "relate", "AB"),

const RELATE_SKIPLIST = Set{Tuple{String, Int, String, String}}([
])
