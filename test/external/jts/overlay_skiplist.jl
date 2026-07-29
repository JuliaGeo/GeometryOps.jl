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
    # "http://trac.osgeo.org/geos/ticket/275" (TestOverlayMisc.xml case 1):
    # STALE EXPECTATION, not a divergence. All three engines agree on the area of
    # the union to 12 significant figures (41227104.2737); the XML expectation is
    # itself 2.96e-12 of area away from what GEOS 3.14 computes today, and ours is
    # 2.15e-10 away — i.e. 5e-18 relative on a 4.1e7 area, pure emission rounding.
    # The positive claim is asserted in xml_suite.jl ("stale XML expectations").
    ("TestOverlayMisc.xml", 1, "union", "AB"),
    # "http://trac.osgeo.org/geos/ticket/488" (TestOverlayMisc.xml case 2):
    # STALE EXPECTATION. This was an engine bug (`found two shells in EdgeRing
    # list`) and is now fixed: the collapsed ring is dropped at emission
    # (`_ring_is_collapsed`). The union is valid and `LG.equals` to GEOS 3.14's
    # exactly, and JTS's own robust corpus agrees with it too
    # (`TestOverlay-geos-488.xml` case 2 states `unionArea = 7.67578758245597E-5`;
    # ours is 7.675787582455967e-5). The XML's expected WKT here is 1.3e-3 of area
    # away from all three, and GEOS fails against it by the same 1.0222e-7 of
    # symmetric difference we do. Asserted in xml_suite.jl.
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
# EVERY case below has `LG.isValidReason(A) == LG.isValidReason(B) == "Valid
# Geometry"` — the corpus is checked for that on the way in (`valid_fn`), and a
# sweep over all 44 files confirmed it for the whole ledger. So none of these is
# an invalid-input case; the engine's §2.2 contract is satisfied throughout.
#
# WHAT THE INPUTS HAVE IN COMMON — sub-ULP slivers in the exact arrangement.
# These are near-coincident-linework regressions: the two inputs carry the same
# real-world boundary digitised twice, a few ULP apart, or a segment passing a
# hair from a vertex. The arrangement is exact, so it *keeps* the resulting
# hair-thin sliver — a corridor or needle whose width is 0–30 ULP of the
# coordinates. JTS does not: its noder computes intersection points as Float64
# while noding, so the two sides of such a sliver land on the same coordinates,
# `EdgeMerger` merges them, their depth deltas cancel, and `Edge.labelDim` labels
# the merged edge `DIM_COLLAPSE` — which `markResultAreaEdges` skips and
# `LineBuilder` emits as a result line. (GEOS's answers show exactly that:
# `intersection` on `TestOverlay-geos-275` case 1 and `TestOverlay-jts-798`
# case 1 is a MULTILINESTRING, not a polygon.)
#
# WHAT THAT COSTS US, NOW — exactly one symptom, and it is a *linework* symptom:
# the emitted result has the right AREA but fails `isValid`. Every entry below is
# `invalid result geometry from …`, and every GEOS reason behind those entries is
# `Self-intersection` or `Ring Self-intersection` at a sliver apex: a corridor
# 0–30 ULP wide rounds to linework that touches itself, or to two rings that
# touch each other, which OGC forbids even though the region it bounds is
# correct. The area identities hold across the whole corpus to 1.4e-11 relative
# (worst: `TestOverlay-qgis-29400` case 4), and no case raises an error.
#
# WHERE THE LEDGER USED TO BE WRONG — this header previously blamed the whole
# cluster on emission-time sliver *representation*, and listed
# `OverlayTopologyError: unable to assign free hole to a shell` / `found two
# shells in EdgeRing list` (51 op-runs, 14 whole cases) plus two genuinely wrong
# areas among its symptoms. Those were NOT representation: they were a design §0
# violation in result assembly. `_compute_ring!` took each minimal ring's
# shell-vs-hole role from `Orientation.isCCW` over the ring's *emitted* Float64
# coordinates, so a ring thinner than an ULP — whose rounded image is inverted or
# self-touching — was filed with the opposite role, and the builder then found
# either two shells in one maximal ring or a hole no shell contains. The role is
# now decided by the sign of the ring's exact signed area over the arrangement's
# node coordinates (`_ring_is_ccw_exact`, maximal_edge_ring.jl). All 51 errors
# are gone, `TestOverlay-misc-1` case 1 and `TestOverlay-misc-4` case 1 now pass
# outright, and the two "area genuinely wrong" entries — `TestOverlay-pg-2055`
# case 1 (`difference` 1395.9967 → 1334.0032858823013, GEOS-identical) and
# `TestOverlay-qgis-37032` case 3 (`union` 1106.2171 → 1104.868307403492 against
# GEOS's 1104.8683074034275) — are cured exactly.
#
# The exactly-representable half of the representation class is fixed too: a ring
# whose emitted coordinates carry repeated points (JTS
# `CoordinateList.add(pt, false)`, which the port was missing) or that reduces to
# fewer than three distinct emitted vertices is dropped from the area result —
# `_ring_is_collapsed` in maximal_edge_ring.jl, the emission-time analogue of
# JTS's `DIM_COLLAPSE`. That flipped `TestOverlay-geos-488` case 1,
# `TestOverlay-jts-808` case 1 and `TestOverlay-misc-4` case 2. What remains
# below is the part where the two sliver sides round to coordinates that are
# *close but not equal* (1–30 ULP), where no exact test can distinguish the
# sliver from a genuine thin polygon. Closing that gap needs a decision about the
# design premise (§0) — it is what JTS's snapping/snap-rounding ladder is for —
# and is deliberately NOT done here.
const ROBUST_KNOWN_DEFECTS = Set{Tuple{String, Int}}([
    ("TestOverlay-geos-1051.xml", 1),                    # invalid result geometry from difference, symdifference
    ("TestOverlay-geos-275.xml", 1),                     # invalid result geometry from symdifference
    ("TestOverlay-geos-350.xml", 1),                     # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-geos-358.xml", 1),                     # invalid result geometry from symdifference
    ("TestOverlay-geos-368.xml", 1),                     # invalid result geometry from union, difference, symdifference
    ("TestOverlay-geos-398.xml", 1),                     # invalid result geometry from symdifference
    ("TestOverlay-geos-522.xml", 1),                     # invalid result geometry from union, difference, symdifference
    ("TestOverlay-geos-615.xml", 1),                     # invalid result geometry from union, symdifference
    ("TestOverlay-geos-737.xml", 1),                     # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-geos-838.xml", 1),                     # invalid result geometry from union, symdifference
    ("TestOverlay-geos-997-union-fail.xml", 1),          # invalid result geometry from difference, symdifference
    ("TestOverlay-geos-list.xml", 1),                    # invalid result geometry from difference, symdifference
    ("TestOverlay-jts-300.xml", 1),                      # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-jts-798.xml", 1),                      # invalid result geometry from union, difference, symdifference
    ("TestOverlay-jts-798.xml", 2),                      # invalid result geometry from union, difference, symdifference
    ("TestOverlay-jts-798.xml", 3),                      # invalid result geometry from union, difference, symdifference
    ("TestOverlay-misc-1.xml", 2),                       # invalid result geometry from symdifference
    ("TestOverlay-misc-1.xml", 3),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-1.xml", 4),                       # invalid result geometry from symdifference
    ("TestOverlay-misc-1.xml", 5),                       # invalid result geometry from symdifference
    ("TestOverlay-misc-2.xml", 1),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-2.xml", 2),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-2.xml", 3),                       # invalid result geometry from intersection, symdifference
    ("TestOverlay-misc-2.xml", 4),                       # invalid result geometry from intersection, symdifference
    ("TestOverlay-misc-2.xml", 5),                       # invalid result geometry from intersection, symdifference
    ("TestOverlay-misc-2.xml", 6),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-2.xml", 7),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-3.xml", 1),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-misc-3.xml", 2),                       # invalid result geometry from union, symdifference
    ("TestOverlay-misc-3.xml", 3),                       # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-misc-3.xml", 4),                       # invalid result geometry from symdifference
    ("TestOverlay-misc-4.xml", 5),                       # invalid result geometry from difference, symdifference
    ("TestOverlay-osmwater.xml", 1),                     # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-osmwater.xml", 2),                     # invalid result geometry from difference, symdifference
    # Was the corpus's worst residual (61.993 on a total of 1396.0, 4.44e-2
    # relative) and is now exact: `difference` is 1334.0032858823013, GEOS's
    # answer to the last bit. A is a rectangle whose "hole" is a needle — its two
    # legs differ by ~4e-6 in one coordinate — and B is a triangle whose apex sits
    # 1.5e-8 from the needle's tip. The B-shaped ring used to come out oriented as
    # a shell rather than a hole of A's rectangle, adding B's 31.0 of area instead
    # of subtracting it; the exact role decision fixes that. What is left is the
    # needle's emitted linework: `Ring Self-intersection` at the needle tip.
    ("TestOverlay-pg-2055.xml", 1),                      # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-pg-2176.xml", 1),                      # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-pg-4182-2.xml", 1),                    # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-pg-4538.xml", 1),                      # invalid result geometry from difference, symdifference
    ("TestOverlay-pg-list.xml", 3),                      # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-pg-list.xml", 5),                      # invalid result geometry from union, intersection, difference, symdifference
    ("TestOverlay-qgis-29400.xml", 1),                   # invalid result geometry from symdifference
    ("TestOverlay-qgis-29400.xml", 2),                   # invalid result geometry from difference, symdifference
    ("TestOverlay-qgis-29400.xml", 3),                   # invalid result geometry from difference, symdifference
    ("TestOverlay-qgis-29400.xml", 4),                   # invalid result geometry from symdifference — and the corpus's worst area identity residual, 1.3e-11 relative
    ("TestOverlay-qgis-29400.xml", 6),                   # invalid result geometry from symdifference
    ("TestOverlay-qgis-37032.xml", 2),                   # invalid result geometry from union, symdifference
    # Same story as pg-2055 case 1: the identity residual was 1.3488 on a total of
    # 1104.87 (1.22e-3 relative) from the same needle-driven role flip, and
    # `union` is now 1104.868307403492 against GEOS's 1104.8683074034275. The
    # remaining failure is `Self-intersection` in the emitted linework.
    ("TestOverlay-qgis-37032.xml", 3),                   # invalid result geometry from union, symdifference
    ("TestOverlay-qgis-37032.xml", 4),                   # invalid result geometry from intersection, difference, symdifference
    ("TestOverlay-qgis-37032.xml", 5),                   # invalid result geometry from difference, symdifference
    ("TestOverlay-rsf-794.xml", 1),                      # invalid result geometry from symdifference
])
