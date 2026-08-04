# Tests for the OverlayNG phase-1 noding substrate (design §2.9): the
# `NodedArrangement` invariants, along-segment ordering vs the exact authority,
# the certified emission fast paths, and the rounded-arrangement / classification
# censuses on a small Natural Earth subset.

using Test
include(joinpath(@__DIR__, "common.jl"))
import GeometryOps: Planar, Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint, UnitSphereFromGeographic
using LinearAlgebra: cross, dot, norm

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_crossing_ids(arr) = [Int32(i) for i in 1:GO.num_nodes(arr) if arr.nodes.keys[i].is_crossing]

# every proper crossing appears as one shared node id on exactly one A-segment and
# one B-segment interior list (invariant 1); node ids are unique (invariant 2);
# no NodedEdge is zero-length (invariant 3).
function check_invariants(arr, na_strings)
    @test length(unique(arr.nodes.keys)) == GO.num_nodes(arr)          # 2
    @test all(e -> e.node_lo != e.node_hi, arr.edges)                  # 3
    @test all(_crossing_ids(arr)) do cid                               # 1
        a_hits = b_hits = 0
        for ((si, _), ids) in arr.seg_nodes
            cid in ids && (si <= na_strings ? (a_hits += 1) : (b_hits += 1))
        end
        a_hits >= 1 && b_hits >= 1
    end
end

# exact-always along-parameter order of a segment's interior node ids
function exact_order(arr, s0, s1, ids)
    R = Rational{BigInt}
    dxr = R(GI.x(s1)) - R(GI.x(s0)); dyr = R(GI.y(s1)) - R(GI.y(s0))
    param(id) = let p = GO._exact_node_point(arr.nodes.keys[id])
        (p[1] - R(GI.x(s0))) * dxr + (p[2] - R(GI.y(s0))) * dyr
    end
    return sort(ids; by = param)
end

# exact-always along-parameter order on the sphere
function exact_order_sph(arr, s0, s1, ids)
    Ne = GO._cross3(GO._vec3(True(), s0), GO._vec3(True(), s1))
    dir(id) = GO._exact_node_dir(True(), arr.nodes.keys[id])
    return sort(ids; lt = (i, j) -> GO._dot3(GO._cross3(dir(i), dir(j)), Ne) > 0)
end

# ---------------------------------------------------------------------------
# 1. Invariants on constructed cases (§2.1)
# ---------------------------------------------------------------------------

@testset "two crossing quads" begin
    A = GI.Polygon([[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 2.0), (6.0, 2.0), (6.0, 6.0), (2.0, 6.0), (2.0, 2.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        na = count(ss -> ss.is_a, arr.segstrings)
        @test length(_crossing_ids(arr)) == 2
        check_invariants(arr, na)
        if m isa Planar
            cpts = sort([GO.node_point(arr, i) for i in _crossing_ids(arr)])
            @test cpts == [(2.0, 4.0), (4.0, 2.0)]     # bit-exact on the integer grid
        end
    end
end

@testset "degree-6 node (tier-2 merge of distinct keys)" begin
    # two A lines and one B line all through the origin -> THREE distinct crossing
    # keys coincident there: crossing_node(lineA1,B), crossing_node(lineA2,B), and
    # crossing_node(lineA1,lineA2) from A's own self-noding. Tier 2 must merge all
    # three into one node — which is what makes the node degree 6.
    A = GI.MultiLineString([[(-1.0, -1.0), (1.0, 1.0)], [(-1.0, 1.0), (1.0, -1.0)]])
    B = GI.LineString([(-1.0, 0.0), (1.0, 0.0)])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        cids = _crossing_ids(arr)
        @test length(cids) == 1                          # merged into one node
        p = GO.node_point(arr, cids[1])
        @test isapprox(p[1], 0.0; atol = 1e-12) && isapprox(p[2], 0.0; atol = 1e-12)
        # the merged node is incident to both A lines and B (three parent strings)
        na = count(ss -> ss.is_a, arr.segstrings)
        check_invariants(arr, na)
    end
end

@testset "crossing exactly on a third string's vertex" begin
    # A horizontal line crosses B-line-1 (vertical) at the origin, which is also
    # B-line-2's endpoint vertex: the crossing key and the vertex key coincide
    # and tier 2 merges them.
    A = GI.LineString([(-2.0, 0.0), (2.0, 0.0)])
    B = GI.MultiLineString([[(0.0, -2.0), (0.0, 2.0)], [(0.0, 0.0), (1.0, 1.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        # the origin is a single node shared by A, B-line-1 and B-line-2's vertex
        origin_ids = [i for i in 1:GO.num_nodes(arr)
            if isapprox(GO.node_point(arr, i)[1], 0.0; atol = 1e-12) &&
               isapprox(GO.node_point(arr, i)[2], 0.0; atol = 1e-12)]
        @test length(origin_ids) == 1
    end
end

@testset "collinear shared boundary — zero phantom crossings" begin
    # edge-adjacent squares sharing the vertex-identical edge x = 2
    A = GI.Polygon([[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]])
    B = GI.Polygon([[(2.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 2.0), (2.0, 0.0)]])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())
        #-- zero phantom crossings: a crossing node is interned only alongside
        #-- its two interior records, so an empty `seg_nodes` implies none exists
        @test sum(length, values(arr.seg_nodes); init = 0) == 0
    end
end

@testset "a-b-a spike input" begin
    # B retraces (0,0)->(1,1)->(0,0); A crosses the retraced segment at one point,
    # reported by both candidate pairs but the same canonical crossing key.
    A = GI.LineString([(-1.0, 0.5), (2.0, 0.5)])
    B = GI.LineString([(0.0, 0.0), (1.0, 1.0), (0.0, 0.0)])
    for m in (Planar(), Spherical())
        arr = GO.NodedArrangement(m, A, B; exact = True())      # must not throw
        cids = _crossing_ids(arr)
        @test length(cids) == 1                                 # one merged node
        m isa Planar && @test GO.node_point(arr, cids[1]) == (0.5, 0.5)
    end
end

# ---------------------------------------------------------------------------
# 2. Ordering cross-check: float-filtered order == exact-always order (§2.5)
# ---------------------------------------------------------------------------

@testset "dense comb ordering matches exact (planar)" begin
    A = GI.LineString([(0.0, 0.0), (201.0, 0.0)])
    B = GI.MultiLineString([[(Float64(i) + 0.3, -1.0), (Float64(i) + 0.3, 1.0)] for i in 1:200])
    arr = GO.NodedArrangement(Planar(), A, B; exact = True())
    # A is string 1, its single segment carries all 200 interior crossings
    ids = arr.seg_nodes[(Int32(1), Int32(1))]
    @test length(ids) == 200
    s0 = arr.segstrings[1].pts[1]; s1 = arr.segstrings[1].pts[2]
    @test ids == exact_order(arr, s0, s1, ids)          # elementwise
    #-- strictly increasing along the segment. The line above already checks the
    #-- order against an independent `Rational{BigInt}` oracle; this adds only
    #-- comparator antisymmetry, so it does not need 199 assertions to say it.
    @test all(c -> GO.rk_compare_along_segment(Planar(), s0, s1,
                  arr.nodes.keys[ids[c-1]], arr.nodes.keys[ids[c]]; exact = True()) < 0,
              2:length(ids))
end

@testset "dense comb ordering matches exact (spherical)" begin
    A = GI.LineString([(0.0, 0.0), (60.0, 0.0)])
    B = GI.MultiLineString([[(Float64(i) * 0.25 + 0.1, -1.0), (Float64(i) * 0.25 + 0.1, 1.0)] for i in 1:200])
    arr = GO.NodedArrangement(Spherical(), A, B; exact = True())
    ids = arr.seg_nodes[(Int32(1), Int32(1))]
    @test length(ids) == 200
    s0 = arr.segstrings[1].pts[1]; s1 = arr.segstrings[1].pts[2]
    @test ids == exact_order_sph(arr, s0, s1, ids)
end

# ---------------------------------------------------------------------------
# 3. Emission certificate audit (§2.6)
# ---------------------------------------------------------------------------

@testset "planar emission: certified == rational, every node" begin
    # dense generic-slope grid + shifted-self coastline-like crossings
    Ag = GI.MultiLineString([[(Float64(k) * 4.0, 0.0), (Float64(k) * 4.0 + 0.31, 1000.0)] for k in 1:60])
    Bg = GI.MultiLineString([[(0.0, Float64(j) * 4.0), (1000.0, Float64(j) * 4.0 + 0.29)] for j in 1:60])
    arr = GO.NodedArrangement(Planar(), Ag, Bg; exact = True())
    #-- accumulated, the way the spherical sibling below already does it: the
    #-- per-node form asserted 7 202 times to establish three properties
    ncert = 0; ntot = 0; cert_ok = true; emit_ok = true
    for i in _crossing_ids(arr)
        k = arr.nodes.keys[i]
        (x, y, cert) = GO._certified_crossing(k.pt, k.a1, k.b0, k.b1)
        rx, ry = GO._exact_crossing_point(k)
        rat = (Float64(rx), Float64(ry))
        ntot += 1
        cert && (ncert += 1; ((x, y) == rat) || (cert_ok = false))
        (GO.node_point(arr, i) == rat) || (emit_ok = false)
    end
    @test ntot > 1000
    @test ncert == ntot                               # 100% certified on clean data (S3)
    @test cert_ok                                     # certified => equal to the rational answer
    @test emit_ok                                     # node_point is the rational answer either way
end

@testset "spherical emission: direction within bound of exact" begin
    Ag = GI.MultiLineString([[(Float64(k) * 0.09 + 0.05, 0.0), (Float64(k) * 0.09 + 0.05 + 0.031, 20.0)] for k in 1:60])
    Bg = GI.MultiLineString([[(0.0, Float64(j) * 0.09 + 0.05), (20.0, Float64(j) * 0.09 + 0.05 + 0.029)] for j in 1:60])
    arr = GO.NodedArrangement(Spherical(), Ag, Bg; exact = True())
    maxdev = 0.0
    for i in _crossing_ids(arr)
        k = arr.nodes.keys[i]
        emitted = GO.node_point(arr, i)
        exact = GO._dir_to_lonlat(GO._sph_crossing_dir(True(), k))
        maxdev = max(maxdev, abs(emitted[1] - exact[1]), abs(emitted[2] - exact[2]))
    end
    @test length(_crossing_ids(arr)) > 1000
    @test maxdev <= 1e-8                              # measured ≤1.4e-14° (S3)
end

# ---------------------------------------------------------------------------
# 4 & 5. Natural Earth subset: rounded-arrangement audit + classification census
# ---------------------------------------------------------------------------

# Rounded-arrangement audit — the one assertion in this file about the property
# the whole substrate exists to guarantee: rounding at emission does not
# introduce topology. No edge incident to a crossing node may properly cross an
# OPPOSITE-SIDE edge once both are realized as Float64.
#
# The scope has to reach PAST the shared node. This function used to compare only
# edges incident to the *same* crossing node, and every such pair carries that
# node's single emitted coordinate as an endpoint — two segments sharing an
# endpoint are `SS_TOUCH` by construction and can never be `SS_PROPER`. Measured
# on the five fixtures below, the old audit examined 184 / 88 / 36 / 136 / 200
# pairs and every one of them shared an endpoint, so it returned 0 by
# construction and would have returned 0 for any arrangement whatsoever.
#
# The corrected audit is discriminating: with emission snapped to 0.1° it finds
# 7 / 0 / 2 / 4 / 8 introduced crossings on the same fixtures, where the old one
# still reports 0.
#
# `m` is the manifold the emitted coordinates are interpreted in, and it is the
# arrangement's own — a spherical arrangement emits lon/lat that downstream
# stages read as spherical, so a lon/lat-planar audit of it would be asking a
# different question than the one the guarantee is about. Emission is lon/lat on
# both manifolds, so the spherical audit lifts back to the kernel's unit vectors
# first; that lift is exactly what the graph builder does with the same points.
emitted_pt(::Planar, p) = p
emitted_pt(::Spherical, p) = UnitSphericalPoint(GI.PointTrait(), p)

function rounded_crossings(m, arr)
    ends = [(emitted_pt(m, GO.node_point(arr, e.node_lo)),
             emitted_pt(m, GO.node_point(arr, e.node_hi))) for e in arr.edges]
    is_a = [arr.segstrings[e.string_idx].is_a for e in arr.edges]
    incident = [n for (n, e) in enumerate(arr.edges)
                if arr.nodes.keys[e.node_lo].is_crossing || arr.nodes.keys[e.node_hi].is_crossing]
    introduced = 0
    for n in incident
        (pa0, pa1) = ends[n]
        for mi in eachindex(arr.edges)
            mi == n && continue
            is_a[mi] == is_a[n] && continue          # opposite side only
            (pb0, pb1) = ends[mi]
            GO.rk_classify_intersection(m, pa0, pa1, pb0, pb1; exact = True()).kind ==
                GO.SS_PROPER && (introduced += 1)
        end
    end
    return introduced
end

ne_ok = try
    import NaturalEarth, GeoJSON
    include(joinpath(@__DIR__, "..", "..", "..", "data", "natural_earth_pairs.jl"))
    global ne_names, ne_geoms = load_ne(110)
    length(ne_geoms) > 0
catch err
    @info "Natural Earth subset skipped (data unavailable)" err
    false
end

# The flag check this block used to run — "every touch/collinear classification
# carries a vertex-incidence flag", via a `classify_census` helper that
# re-enumerated the A×B candidate pairs with the same index, the same traversal
# and the same kernel call the noder uses — is a live `@assert` inside
# `collect.jl` three lines from the classification itself. A test-side
# reimplementation of noding stage 1 is not a second opinion; it is the same
# opinion, computed twice.
@testset "Natural Earth subset (rounded-arrangement audit)" begin
    if !ne_ok
        @test_skip "Natural Earth data unavailable"
    else
        picks = String["Brazil", "France", "Egypt", "India", "Australia"]
        tested = 0
        for nm in picks
            idx = findfirst(==(nm), ne_names)
            idx === nothing && continue
            A = ne_geoms[idx]
            B = shift_geom(A, 0.5, 0.0)
            tested += 1
            arr_p = GO.NodedArrangement(Planar(), A, B; exact = True())
            arr_s = GO.NodedArrangement(Spherical(), A, B; exact = True())
            #-- the fixture must actually cross, or the audit below is vacuous
            @test !isempty(_crossing_ids(arr_p)) && !isempty(_crossing_ids(arr_s))
            @test rounded_crossings(Planar(), arr_p) == 0
            @test rounded_crossings(Spherical(), arr_s) == 0
        end
        @test tested >= 2
    end
end
