# Skiplist and known-defect ledger for the JTS overlay XML suites
# (`overlay_runner.jl`, driven from `test/methods/clipping/overlayng/xml_suite.jl`).
#
# `OVERLAY_SKIPLIST` entries are `(file, case_index, op name, arg order)` tuples
# identifying one `<op>` in one `<case>` of one vendored XML file, in exactly the
# shape used by `relate_skiplist.jl`:
#
# - `file` is the basename, e.g. "TestOverlayAA.xml".
# - `case_index` is the 1-based index of the `<case>` within the file.
# - `op name` is as written in the XML, e.g. "differenceNG".
# - `arg order` is "AB" when the op's `arg1` is the case's A geometry and "BA"
#   when it is B (the `*NG` files run `difference` both ways round).
#
# DISCIPLINE: every entry MUST carry a comment giving the case description and
# explaining exactly why GeometryOps diverges from the XML expectation (or why
# the case cannot run). Entries without a justification must not be merged.
# Skipped ops are reported by the runner — they are never silently dropped.
#
# Not listed here, because the runner skips them structurally with a recorded
# reason: point inputs (phase 3), GEOMETRYCOLLECTION inputs, ops that are not one
# of the four overlay ops, and runs with a FIXED precision model.

const OVERLAY_SKIPLIST = Set{Tuple{String, Int, String, String}}([
    # "mLmA - A and B complex, overlapping and touching" (TestOverlayLA.xml case
    # 2): ENGINE BUG. A is a donut whose hole contains a second donut; B is a
    # MULTILINESTRING one of whose components runs along the outer donut's hole
    # boundary (y = 240, x ∈ [80, 200]) while the other crosses it at x = 120 and
    # x = 140. The three sub-segments of that hole-boundary-collinear line are
    # dropped from the intersection (our result has length 500, JTS's and GEOS's
    # 620). With either B component alone the answer is correct, so the trigger
    # is a line collinear with a hole boundary that is *split* by nodes coming
    # from another input component. Reproducer in the report; smallest form is
    # A = the first MULTIPOLYGON member alone with the full B.
    ("TestOverlayLA.xml", 2, "intersection", "AB"),
    # "http://trac.osgeo.org/geos/ticket/275" (TestOverlayMisc.xml case 1):
    # STALE EXPECTATION, not a divergence. All three engines agree on the area of
    # the union to 12 significant figures (41227104.2737); the XML expectation is
    # itself 2.96e-12 of area away from what GEOS 3.14 computes today, and ours is
    # 2.15e-10 away — i.e. 5e-18 relative on a 4.1e7 area, pure emission rounding.
    # The positive claim is asserted in xml_suite.jl ("stale XML expectations").
    ("TestOverlayMisc.xml", 1, "union", "AB"),
    # "http://trac.osgeo.org/geos/ticket/488" (TestOverlayMisc.xml case 2):
    # ENGINE BUG — `OverlayTopologyError: found two shells in EdgeRing list` on
    # the union of two VALID multipolygons. Same collapsed-result-ring class as
    # the robust ledger below (see ROBUST_KNOWN_DEFECTS).
    ("TestOverlayMisc.xml", 2, "union", "AB"),
    # "https://trac.osgeo.org/geos/ticket/368" (TestOverlayMisc.xml case 3):
    # STALE EXPECTATION. The XML's expected WKT is written to 8 decimal places
    # (e.g. `-199983.26344477` where the true coordinate is
    # `-199983.26344477304`), which alone accounts for the 3.03e-7 area
    # difference. Our result is `equals` to GEOS 3.14's exactly (symmetric
    # difference area 0.0), and GEOS also fails against this expectation.
    ("TestOverlayMisc.xml", 3, "intersection", "AB"),
    # "https://trac.osgeo.org/geos/ticket/522" (TestOverlayMisc.xml case 4):
    # EMISSION ROUNDING, ours is at least as good. The two results differ only in
    # the last 2 ulps of one y coordinate (851253.4627870634 vs ...636 at 8.5e5,
    # ulp = 1.16e-10). Ours is 1.37e-9 of area from the XML expectation, GEOS is
    # 5.93e-9 from it — i.e. we are closer to the stated answer than GEOS is.
    ("TestOverlayMisc.xml", 4, "intersection", "AB"),
    # "https://trac.osgeo.org/geos/ticket/737" (TestOverlayMisc.xml case 5):
    # INVALID INPUT. `LG.isValidReason(A)` is
    # `Self-intersection[374476.386125697 4787775.46792226]`; the engine
    # contracts on valid input (design §2.2) and raises an OverlayTopologyError.
    ("TestOverlayMisc.xml", 5, "union", "AB"),
])

# ## Known-defect ledger for the `robust/overlay` corpus
#
# The `test/data/jts/overlay_robust/` corpus is JTS's harvest of overlay
# regressions from GEOS / PostGIS / QGIS / shapely / JTS bug reports: hard, valid,
# real-world inputs. `xml_suite.jl` runs the oracle-free area-identity sweep over
# all of it and pins the exact set of `(file, case_index)` pairs that currently
# fail, so both a regression and a fix show up as a set difference.
#
# ONE root cause dominates the ledger (confirmed on the 8-point reproducer
# `TestOverlay-jts-798.xml` case 1, quoted in full in `xml_suite.jl`):
#
#   COLLAPSED RESULT RINGS. When part of the result collapses to zero width — a
#   sliver whose two sides round to the same emitted segment, i.e. an
#   out-and-back spike — the engine keeps it inside an area ring instead of
#   excluding it from the area result and emitting it as a result *line* the way
#   JTS's OverlayEdgeRing/LineBuilder split does. The visible symptoms are
#     * `isValid` false with GEOS reason `Too few points in geometry component`,
#       `Self-intersection` or `Ring Self-intersection` at the spike apex;
#     * `OverlayTopologyError: unable to assign free hole to a shell` when the
#       *whole* result collapses (there are then holes but no shell at all);
#     * `OverlayTopologyError: found two shells in EdgeRing list` when a collapsed
#       ring is classified as a second shell of one maximal ring.
#   In every ledger case where an area *is* produced it agrees with GEOS to the
#   rounding band, so this is a result-*representation* defect, not an arithmetic
#   one. The two exceptions, where the area is genuinely wrong, are called out
#   individually below.
const ROBUST_KNOWN_DEFECTS = Set{Tuple{String, Int}}([
    ("TestOverlay-geos-1051.xml", 1),                    # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-geos-275.xml", 1),                     # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-geos-350.xml", 1),                     # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-geos-358.xml", 1),                     # invalid result geometry from union, symdifference
    ("TestOverlay-geos-368.xml", 1),                     # invalid result geometry from union, difference, symdifference
    ("TestOverlay-geos-398.xml", 1),                     # OverlayTopologyError: found two shells in EdgeRing list
    ("TestOverlay-geos-488.xml", 1),                     # OverlayTopologyError: found two shells in EdgeRing list
    ("TestOverlay-geos-522.xml", 1),                     # invalid result geometry from union, difference, symdifference
    ("TestOverlay-geos-615.xml", 1),                     # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-geos-737.xml", 1),                     # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-geos-838.xml", 1),                     # invalid result geometry from union, intersection, symdifference
    ("TestOverlay-geos-997-union-fail.xml", 1),          # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-geos-list.xml", 1),                    # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-jts-300.xml", 1),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-jts-798.xml", 1),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-jts-798.xml", 2),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-jts-798.xml", 3),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-jts-808.xml", 1),                      # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-1.xml", 1),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-1.xml", 2),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-1.xml", 3),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-1.xml", 4),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-1.xml", 5),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-2.xml", 1),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-2.xml", 2),                       # invalid result geometry from union, intersection, symdifference
    ("TestOverlay-misc-2.xml", 3),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-2.xml", 4),                       # invalid result geometry from intersection, symdifference
    ("TestOverlay-misc-2.xml", 5),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-2.xml", 6),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-2.xml", 7),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-3.xml", 1),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-3.xml", 2),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-3.xml", 3),                       # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-misc-3.xml", 4),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-4.xml", 1),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-misc-4.xml", 2),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-4.xml", 5),                       # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-osmwater.xml", 1),                     # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-osmwater.xml", 2),                     # invalid result geometry from difference, symdifference
    # AREA WRONG (not just representation): identity residual 61.99 on a total of
    # 1396.0, i.e. 4.44e-2 relative — the worst in the corpus. A is a rectangle
    # whose "hole" is itself a zero-area out-and-back spike (its two legs differ
    # by ~4e-7 in one coordinate); `difference` then returns 1395.997 where GEOS
    # returns 1334.003, double-counting the spike region.
    ("TestOverlay-pg-2055.xml", 1),
    ("TestOverlay-pg-2176.xml", 1),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-pg-4182-2.xml", 1),                    # OverlayTopologyError: found two shells in EdgeRing list
    ("TestOverlay-pg-4538.xml", 1),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-pg-list.xml", 3),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-pg-list.xml", 5),                      # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-qgis-29400.xml", 1),                   # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-qgis-29400.xml", 2),                   # invalid result geometry from difference, symdifference
    ("TestOverlay-qgis-29400.xml", 3),                   # invalid result geometry from difference, symdifference
    ("TestOverlay-qgis-29400.xml", 4),                   # invalid result geometry from intersection, symdifference
    ("TestOverlay-qgis-29400.xml", 6),                   # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-qgis-37032.xml", 2),                   # invalid result geometry from union, symdifference
    # AREA WRONG (not just representation): identity residual 1.349 on a total of
    # 1104.9, i.e. 1.22e-3 relative — same spike-hole shape as pg-2055 case 1.
    ("TestOverlay-qgis-37032.xml", 3),
    ("TestOverlay-qgis-37032.xml", 4),                   # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-qgis-37032.xml", 5),                   # OverlayTopologyError: unable to assign free hole to a shell
    ("TestOverlay-rsf-794.xml", 1),                      # OverlayTopologyError: unable to assign free hole to a shell
])
