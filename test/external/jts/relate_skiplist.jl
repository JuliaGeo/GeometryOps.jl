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
    # "P/L-2: a point and a zero-length line" (validate/TestRelatePL.xml,
    # case 2): JTS RelateNG (and GEOS >= 3.13) deliberately treats zero-length
    # LineStrings as topologically identical to Points (RelateNG.java class
    # javadoc, difference 3), so `equalsTopo(POINT(110 200),
    # LINESTRING(110 200, 110 200))` is `true`. The legacy XML expectation
    # (`false`) encodes old-RelateOp behavior, where the line's declared
    # dimension (1) differs from the point's (0). GEOS 3.14.1 RelateNG agrees
    # with GeometryOps. (general/TestRelatePL.xml shares this basename, but its
    # case 2 has only `relate` ops, so the key is unambiguous.)
    ("TestRelatePL.xml", 2, "equalsTopo", "AB"),
    # "A/P - Point is on boundary of polygon" (validate/TestRobustRelateFloat.xml,
    # case 1): the case intends POINT(0.95 0.05) to lie on the hypotenuse of
    # POLYGON((0 0, 1 0, 0 1, 0 0)), but in exact double semantics
    # 0.95 + 0.05 == 1 - 3*2^-56 < 1, so the parsed point lies strictly inside
    # the polygon and `contains` is exactly `true`. The legacy expectation
    # (`false`) assumes the idealized, non-representable coordinates. GEOS
    # 3.14.1 RelateNG also returns `true`.
    ("TestRobustRelateFloat.xml", 1, "contains", "AB"),
])
