# Skiplist for the JTS relate XML conformance suite (`relate_runner.jl`).
#
# Each entry is a `(file, case description, op name)` tuple identifying one
# `<op>` in one `<case>` of one vendored XML file (file is the basename, e.g.
# "TestRelateAA.xml"; op name is as written in the XML, e.g. "equalsTopo").
#
# DISCIPLINE: every entry MUST be accompanied by a comment explaining exactly
# why GeometryOps diverges from JTS on that case (or why the case cannot run),
# ideally with an issue/task reference. Entries without a justification comment
# must not be merged. Skipped cases are reported by `run_relate_cases` — they
# are never silently dropped.
#
# Example entry:
#     # GO returns FF* for empty-geometry boundary, JTS expects F0*; tracked in #XXXX
#     ("TestRelateEmpty.xml", "P/L - empty point VS empty line", "relate"),

const RELATE_SKIPLIST = Set{Tuple{String, String, String}}([
])
