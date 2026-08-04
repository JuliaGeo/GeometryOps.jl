# Shared fixtures and helpers for the OverlayNG test directory.
#
# Ten test files sweep the same engine, and before this file existed they each
# carried their own copy of the operation tables, the LibGEOS conversion, and
# the 2x2 square pair. The copies had drifted — three different spellings of
# "our result as a LibGEOS geometry", two incompatible `OPS` tuples — so a fix
# to one never reached the others.
#
# Everything here is either a *helper* or one of the two canonical squares.
# Nothing test-specific lives here on purpose: a fixture in this file changes
# the meaning of every file that includes it, so the bar for adding one is that
# it is genuinely shared vocabulary and not just a shape two files happen to
# have in common.
#
# Each file pulls this in with
#
#     include(joinpath(@__DIR__, "common.jl"))
#
# which keeps every file bare-includable and gives each `@safetestset` module
# its own copy of the bindings — the same pattern
# `test/data/natural_earth_pairs.jl` already uses.

import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

# Exact predicates on. Every suite in this directory runs the engine in its
# certified configuration; `exact = False()` is probed in exactly one place
# (`api.jl`'s plumbing test) and nowhere else.
const EX = GO.True()

# ## The canonical pair
#
# Two axis-aligned squares overlapping in the unit square (1,1)-(2,2). Every
# operation on the pair has a hand-checkable answer, which is why this shape is
# the default fixture throughout the directory:
#
#     area(A) = area(B) = 4    A ∩ B = 1    A ∪ B = 7    A ∖ B = 3    A △ B = 6
const SQ_A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
const SQ_B = GI.Polygon([[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0), (1.0, 1.0)]])

# ## Operations
#
# Two spellings are in use: the engine takes `_OverlayOpCode` values, while the
# corpus-driven suites (`xml_suite`, `fuzz`, `s2_differential`) key their tables
# on symbols because that is what the corpora carry.
const OP_CODES = (GO.OVERLAY_INTERSECTION, GO.OVERLAY_UNION,
                  GO.OVERLAY_DIFFERENCE, GO.OVERLAY_SYMDIFFERENCE)
const OP_SYMS = (:intersection, :union, :difference, :symdifference)
const OPCODE = Dict(zip(OP_SYMS, OP_CODES))

opname(op::Symbol) = String(op)
opname(op) = String(OP_SYMS[findfirst(==(op), OP_CODES)])

# ## LibGEOS interop
#
# `to_lg` is idempotent on LibGEOS geometries, so the GEOS helpers below accept
# either side of the boundary and callers need not remember which they hold.
to_lg(g) = g isa LG.AbstractGeometry ? g : GI.convert(LG, g)
const lgc = to_lg          # historical spelling, kept because it reads better inline
wkt(s) = LG.readgeom(s)
giwkt(s) = GO.tuples(LG.readgeom(s))

function geos_op(op, a, b)
    la, lb = to_lg(a), to_lg(b)
    sym = op isa Symbol ? op : Symbol(opname(op))
    sym === :intersection && return LG.intersection(la, lb)
    sym === :union && return LG.union(la, lb)
    sym === :difference && return LG.difference(la, lb)
    return LG.symmetricDifference(la, lb)
end

# Our result → LibGEOS. Empty results become one canonical empty geometry:
# GEOS `equals` is true between any two empties, and the *type* of an empty
# overlay result is what `xml_suite.jl`'s TestOverlayEmpty run checks.
function result_to_lg(g)
    t = GI.trait(g)
    t isa GI.PointTrait && return LG.Point(Float64(GI.x(g)), Float64(GI.y(g)))
    t isa GI.GeometryCollectionTrait &&
        return LG.GeometryCollection(LG.Geometry[result_to_lg(s) for s in GI.getgeom(g)])
    GI.npoint(g) == 0 && return LG.readgeom("GEOMETRYCOLLECTION EMPTY")
    return to_lg(g)
end

# `A` translated by `(dx, dy)`. Overlaying a geometry with a shifted copy of
# itself makes every segment cross its own image, which is the densest crossing
# workload a fixture can have for its vertex count.
shift_geom(g, dx, dy) = GO.apply(GI.PointTrait(), g) do p
    (GI.x(p) + dx, GI.y(p) + dy)
end

# ## The per-dimension signature
#
#     dim_signature(m, g) -> (areal area, lineal length, point-component count)
#
# Every scalar-area oracle in this directory is blind to content of lower
# dimension: a result that silently drops its 1-D components has exactly the
# area of one that keeps them, and a dropped 0-D component moves neither area
# nor length. That is one defect class, not three, and comparing signatures
# componentwise is what closes it — see the mutation `_build_lines -> []`, which
# every area-valued harness scored green.
#
# Polygon boundaries deliberately do NOT contribute to the length slot: only
# lineal components standing on their own do, because a polygon's perimeter is
# a function of its area part and would swamp the very thing being measured.
# (`GO.perimeter` targets `AbstractCurveTrait`, so it walks polygon rings too —
# hence the explicit trait dispatch rather than one call at the top.)
dim_signature(m, g) = _dim_sig(m, GI.trait(g), g)

_dim_sig(m, ::GI.PointTrait, g) = (0.0, 0.0, 1)
_dim_sig(m, ::GI.MultiPointTrait, g) = (0.0, 0.0, GI.ngeom(g))
_dim_sig(m, ::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, g) =
    (Float64(GO.area(m, g)), 0.0, 0)
#-- the `< 2` guard is for empty lineal results: `GO.perimeter` asserts it has
#-- an edge to measure, and an empty component contributes nothing anyway
_dim_sig(m, ::Union{GI.LineStringTrait, GI.LinearRingTrait, GI.MultiLineStringTrait}, g) =
    GI.npoint(g) < 2 ? (0.0, 0.0, 0) : (0.0, Float64(GO.perimeter(m, g)), 0)
function _dim_sig(m, ::GI.GeometryCollectionTrait, g)
    acc = (0.0, 0.0, 0)
    for sub in GI.getgeom(g)
        s = dim_signature(m, sub)
        acc = (acc[1] + s[1], acc[2] + s[2], acc[3] + s[3])
    end
    return acc
end
_dim_sig(m, ::Nothing, g) = (0.0, 0.0, 0)

"""
    signatures_agree(m, ours, theirs; atol, rtol) -> Bool

Compare two results by per-dimension content. The area and length slots are
compared to a tolerance (both engines round their emitted vertices), the point
count exactly (a component is either there or it is not).
"""
function signatures_agree(m, ours, theirs; atol = 1e-9, rtol = 1e-9)
    a, b = dim_signature(m, ours), dim_signature(m, theirs)
    return isapprox(a[1], b[1]; atol, rtol) && isapprox(a[2], b[2]; atol, rtol) && a[3] == b[3]
end
