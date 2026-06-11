# Differential fuzzing of RelateNG against LibGEOS (Task 25).
#
# Property: `string(GO.relate(GO.RelateNG(), a, b)) == GEOSRelate(a, b)` for
# seeded random and adversarial geometry pairs. GEOS >= 3.13 runs RelateNG
# natively, so the oracle is the same algorithm family and any divergence is
# either a real bug on our side or an exactness win (GEOS computes in floating
# point; our kernel is exact via AdaptivePredicates).
#
# Valid-input fuzzing only for polygon inputs: JTS/GEOS RelateNG assumes valid
# polygon topology, so random polygons are filtered through `LG.isValid`.
# Self-intersecting (bowtie) rings are still fuzzed — deliberately, as
# LINESTRING inputs (see `bowtie_pair`), where they are valid geometry.
#
# Case count: ~`FUZZ_N` total, split evenly across the generators below.
# Override with `ENV["GO_RELATENG_FUZZ_N"]` for deep local runs (e.g. 5000).
# Each generator has its own fixed-seed RNG stream, so the first cases of a
# deep run are exactly the CI cases and fixtures stay stable across N.

using Test
using Random
import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG

include(joinpath(@__DIR__, "..", "..", "data", "polygon_generation.jl"))

# The oracle must be RelateNG-native GEOS (3.13+); bail out loudly otherwise.
@assert VersionNumber(LG.GEOS_VERSION) >= v"3.13" "Differential fuzz requires GEOS >= 3.13 (RelateNG-native oracle); found $(LG.GEOS_VERSION)"

const FUZZ_N = something(tryparse(Int, get(ENV, "GO_RELATENG_FUZZ_N", "")), 500)
const ALG = GO.RelateNG()

# -- LibGEOS oracle plumbing --------------------------------------------------

const LG_CTX = LG.get_global_context()
# Full-precision WKT so any printed divergence reproduces the exact doubles
# (the default writer trims to ~16 significant digits, which would erase
# ulp-level perturbations).
const WKT_WRITER = LG.WKTWriter(LG_CTX; trim = true, roundingprecision = 17)

to_lg(g) = GI.convert(LG, g)
full_wkt(lg_geom) = LG.writegeom(lg_geom, WKT_WRITER, LG_CTX)
# LibGEOS has no high-level `relate`; the generated wrapper returns the DE-9IM
# string (and frees the GEOS-owned buffer) itself.
lg_relate(la, lb) = LG.GEOSRelate_r(LG_CTX, la, lb)

# -- Exactness-wins fixtures --------------------------------------------------

# Divergences where our exact kernel is right and FP GEOS is wrong, verified
# by hand with rational arithmetic. Keyed by `(full_wkt(a), full_wkt(b))`;
# the value is OUR (verified-correct) DE-9IM string, asserted instead of the
# oracle. Document the hand verification in a comment on each entry.
const EXACTNESS_WINS = Dict{Tuple{String, String}, String}(
)

# -- The property -------------------------------------------------------------

function check_pair(a, b)
    la, lb = to_lg(a), to_lg(b)
    wa, wb = full_wkt(la), full_wkt(lb)
    ours = try
        string(GO.relate(ALG, a, b))
    catch
        println("RelateNG error on:\n  A = $wa\n  B = $wb")
        rethrow()
    end
    fixture = get(EXACTNESS_WINS, (wa, wb), nothing)
    if fixture !== nothing
        @test ours == fixture
        return nothing
    end
    theirs = lg_relate(la, lb)
    if ours != theirs
        println("RelateNG/LibGEOS divergence:\n  A = $wa\n  B = $wb\n  GO = $ours  GEOS = $theirs")
    end
    @test ours == theirs
    return nothing
end

# -- Random generators (floating-point coordinates) ---------------------------

function random_poly(rng)
    x, y = 4rand(rng) - 2, 4rand(rng) - 2
    nverts = rand(rng, 4:12)
    coords = generate_random_poly(x, y, nverts, 1.0 + rand(rng), 0.3 * rand(rng), 0.2 * rand(rng), rng)
    return GI.Polygon([[(Float64(p[1]), Float64(p[2])) for p in coords[1]]])
end

# `generate_random_poly` does not guarantee non-self-intersection; reject
# invalid candidates (RelateNG assumes valid polygons, see header note).
function random_valid_polygon(rng)
    for _ in 1:50
        poly = random_poly(rng)
        LG.isValid(to_lg(poly)) && return poly
    end
    # Vanishingly unlikely fallback: a unit right triangle at a random offset.
    x, y = 4rand(rng) - 2, 4rand(rng) - 2
    return GI.Polygon([[(x, y), (x + 1, y), (x, y + 1), (x, y)]])
end

function random_linestring(rng; npts = rand(rng, 2:8))
    x, y = 4rand(rng) - 2, 4rand(rng) - 2
    pts = [(x, y)]
    for _ in 2:npts
        x += randn(rng)
        y += randn(rng)
        push!(pts, (x, y))
    end
    return GI.LineString(pts)
end

random_point(rng) = GI.Point(4rand(rng) - 2, 4rand(rng) - 2)

# -- Integer-grid helpers (exact degeneracies: shared vertices, collinear
#    edges, touches actually hit the zero-orientation predicate paths) --------

int_pt(rng) = (Float64(rand(rng, -5:5)), Float64(rand(rng, -5:5)))

# CCW integer triangle ring, optionally through fixed (shared) vertices.
function int_triangle_ring(rng, fixed::Vector{NTuple{2, Float64}} = NTuple{2, Float64}[])
    for _ in 1:100
        pts = copy(fixed)
        while length(pts) < 3
            p = int_pt(rng)
            p in pts || push!(pts, p)
        end
        a, b, c = pts
        cr = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
        cr == 0 && continue
        return cr > 0 ? [a, b, c, a] : [a, c, b, a]
    end
    error("could not build a non-degenerate integer triangle")
end

int_triangle(rng, fixed = NTuple{2, Float64}[]) = GI.Polygon([int_triangle_ring(rng, fixed)])

function int_rect(x0, y0, w, h)
    x0, y0, x1, y1 = Float64(x0), Float64(y0), Float64(x0 + w), Float64(y0 + h)
    return GI.Polygon([[(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]])
end

function diamond(cx, cy, r)
    cx, cy, r = Float64(cx), Float64(cy), Float64(r)
    return GI.Polygon([[(cx - r, cy), (cx, cy - r), (cx + r, cy), (cx, cy + r), (cx - r, cy)]])
end

function int_linestring(rng; npts = rand(rng, 2:6))
    pts = NTuple{2, Float64}[]
    while length(pts) < npts
        p = int_pt(rng)
        (isempty(pts) || pts[end] != p) && push!(pts, p)
    end
    return GI.LineString(pts)
end

# -- Mixed-type random pairs --------------------------------------------------

function line_line_pair(rng)
    return rand(rng, Bool) ? (random_linestring(rng), random_linestring(rng)) :
        (int_linestring(rng), int_linestring(rng))
end

function point_polygon_pair(rng)
    # Integer half: the point frequently lands exactly on a vertex or edge.
    return rand(rng, Bool) ? (random_point(rng), random_valid_polygon(rng)) :
        (GI.Point(int_pt(rng)), int_triangle(rng))
end

function line_polygon_pair(rng)
    return rand(rng, Bool) ? (random_linestring(rng), random_valid_polygon(rng)) :
        (int_linestring(rng), int_triangle(rng))
end

# -- Adversarial generators ---------------------------------------------------

# B reuses 1-2 of A's vertices: exact shared nodes between the inputs.
function shared_vertex_pair(rng)
    ra = int_triangle_ring(rng)
    nshared = rand(rng, 1:2)
    rb = int_triangle_ring(rng, ra[1:nshared])
    return GI.Polygon([ra]), GI.Polygon([rb])
end

# Rectangles sharing the line y = y0: collinear overlapping/abutting edges
# (B below A touches edge-to-edge; B on the same side overlaps interiors).
function collinear_edge_pair(rng)
    x0, y0 = rand(rng, -5:2), rand(rng, -5:2)
    w, h = rand(rng, 1:4), rand(rng, 1:4)
    a = int_rect(x0, y0, w, h)
    dx = rand(rng, -2:w)
    w2, h2 = rand(rng, 1:4), rand(rng, 1:4)
    b = rand(rng, Bool) ? int_rect(x0 + dx, y0 - h2, w2, h2) : int_rect(x0 + dx, y0, w2, h2)
    return a, b
end

# nextfloat/prevfloat by 1-3 ulps, random direction.
nudge(rng, v) = (rand(rng, Bool) ? nextfloat : prevfloat)(v, rand(rng, 1:3))

# Segment configurations sitting within a few ulps of a crossing, touch, or
# collinear overlap — the FP-noise regime where exact and inexact kernels can
# disagree about orientation signs.
function ulp_pair(rng)
    a1 = (2rand(rng) - 1, 2rand(rng) - 1)
    a2 = (a1[1] + (rand(rng) + 0.1) * rand(rng, (-1, 1)), a1[2] + (rand(rng) + 0.1) * rand(rng, (-1, 1)))
    # FP point approximately on A at parameter t (generally a few ulps off).
    seg(t) = (a1[1] + t * (a2[1] - a1[1]), a1[2] + t * (a2[2] - a1[2]))
    maybe_nudge(p) = rand(rng, Bool) ? (nudge(rng, p[1]), rand(rng, Bool) ? nudge(rng, p[2]) : p[2]) : p
    mode = rand(rng, 1:3)
    if mode == 1
        # B ends within ulps of A's interior (near endpoint-on-segment touch).
        m = maybe_nudge(seg(0.2 + 0.6 * rand(rng)))
        b2 = (m[1] + 2rand(rng) - 1, m[2] + 2rand(rng) - 1)
        return GI.LineString([a1, a2]), GI.LineString([m, b2])
    elseif mode == 2
        # B nearly collinear with A, overlapping in parameter range.
        t1, t2 = minmax(2rand(rng) - 0.5, 2rand(rng) - 0.5)
        t1 == t2 && (t2 += 0.5)
        return GI.LineString([a1, a2]), GI.LineString([maybe_nudge(seg(t1)), maybe_nudge(seg(t2))])
    else
        # X-crossing through an ulp-nudged interior point of A.
        m = seg(0.2 + 0.6 * rand(rng))
        d = (-(a2[2] - a1[2]), a2[1] - a1[1])
        s = rand(rng) + 0.1
        b1 = (m[1] - s * d[1], m[2] - s * d[2])
        b2m = maybe_nudge(m)
        b2 = (b2m[1] + s * d[1], b2m[2] + s * d[2])
        return GI.LineString([a1, a2]), GI.LineString([b1, b2])
    end
end

# Zero-length (degenerate) linestrings against everything.
function zero_length_pair(rng)
    mode = rand(rng, 1:4)
    if mode == 1
        p = (2rand(rng) - 1, 2rand(rng) - 1)
        return GI.LineString([p, p]), random_valid_polygon(rng)
    elseif mode == 2
        p = (2rand(rng) - 1, 2rand(rng) - 1)
        return GI.LineString([p, p]), random_linestring(rng)
    elseif mode == 3
        # Zero-length line exactly at a polygon vertex.
        p = int_pt(rng)
        return GI.LineString([p, p]), int_triangle(rng, [p])
    else
        # Two coincident zero-length lines.
        p = int_pt(rng)
        return GI.LineString([p, p]), GI.LineString([p, p])
    end
end

# Polygon rings touching at exactly one point: vertex-vertex (two diamonds)
# or vertex on the interior of an edge (diamond tip on a rectangle side).
function touching_rings_pair(rng)
    if rand(rng, Bool)
        cx, cy = rand(rng, -3:3), rand(rng, -3:3)
        r1, r2 = rand(rng, 1:3), rand(rng, 1:3)
        return diamond(cx, cy, r1), diamond(cx + r1 + r2, cy, r2)
    else
        x0, y0 = rand(rng, -4:1), rand(rng, -4:0)
        w, h = rand(rng, 1:3), rand(rng, 2:4)
        ty = rand(rng, (y0 + 1):(y0 + h - 1))
        r = rand(rng, 1:3)
        return int_rect(x0, y0, w, h), diamond(x0 + w + r, ty, r)
    end
end

# Self-intersecting (bowtie) rings are invalid polygons, so they are fuzzed
# as closed LINESTRINGs, which exercises the line self-noding paths.
function bowtie_ring(rng)
    x0, y0 = Float64(rand(rng, -4:2)), Float64(rand(rng, -4:2))
    s = Float64(rand(rng, 1:3))
    return [(x0, y0), (x0 + 2s, y0 + 2s), (x0 + 2s, y0), (x0, y0 + 2s), (x0, y0)]
end

function bowtie_pair(rng)
    a = GI.LineString(bowtie_ring(rng))
    b = rand(rng, Bool) ? GI.LineString(bowtie_ring(rng)) : random_valid_polygon(rng)
    return a, b
end

# -- Driver -------------------------------------------------------------------

const GENERATORS = [
    ("polygon/polygon", rng -> (random_valid_polygon(rng), random_valid_polygon(rng))),
    ("line/line", line_line_pair),
    ("point/polygon", point_polygon_pair),
    ("line/polygon", line_polygon_pair),
    ("shared vertices", shared_vertex_pair),
    ("collinear edges", collinear_edge_pair),
    ("ulp near-crossings", ulp_pair),
    ("zero-length lines", zero_length_pair),
    ("rings touching at points", touching_rings_pair),
    ("bowtie linestrings", bowtie_pair),
]

const CASES_PER_GENERATOR = max(1, cld(FUZZ_N, length(GENERATORS)))

@testset "RelateNG vs LibGEOS fuzz (N=$FUZZ_N, $CASES_PER_GENERATOR per generator)" begin
    for (i, (name, gen)) in enumerate(GENERATORS)
        rng = Xoshiro(0x9E1A7E20260000 + i)
        @testset "$name" begin
            for _ in 1:CASES_PER_GENERATOR
                a, b = gen(rng)
                check_pair(a, b)
            end
        end
    end
end
