# # Spherical leaves of the lightweight geometry-relation processors

#=
The processors in `geom_geom_processors.jl` carry the DE-9IM allow/require
bookkeeping for `intersects`, `disjoint`, `within`, `contains`, `covers`,
`coveredby` and `touches`. That bookkeeping is manifold-independent; only four
geometric leaves are not, and this file supplies their `Spherical` methods:

| leaf | question |
|:--|:--|
| `_point_segment_orientation` | is a point on a segment, and is it an endpoint? |
| `_point_filled_curve_orientation` | is a point in, on, or out of a filled ring? |
| `_seg_seg_orientation` | how do two segments meet? |
| `_split_segment_interactions` | does a hinging segment run inside, outside, or both? |

Everything here is allocation-free per call. Two properties buy that. Vertices
are converted to `UnitSphericalPoint` lazily through
[`SphericalRingPoints`](@ref) rather than materialized into a `Vector`, so a
predicate against an unprepared ring still touches no heap. And the segment
classifier is symbolic — it names the kind of meeting and which vertices are
incident, and never constructs an intersection coordinate. Constructing one is
what makes great-circle segment intersection awkward (two antipodal candidate
points, one of them wrong); not needing one removes the problem from the
predicate path entirely. The one place a coordinate *is* needed — splitting a
hinging segment — builds it only after the crossing is known to be proper, which
is what makes choosing between the two candidates safe there.

## Which tolerance regime, and why not the exact kernel's

Every predicate here is the tolerance-banded `UnitSpherical` one:
`spherical_orient`, whose `16 * eps` relative band is ~3.5e-15 radians, over a
span test (`_sph_on_arc`) that is correct for arcs of any length. That is deliberate, and it is not the same choice RelateNG
makes.

RelateNG's `rk_*` kernel is exact-or-nothing: at `exact = True()` it decides
degeneracies over `Rational{BigInt}`, and at `exact = False()` it degrades to a
raw Float64 triple product still compared `== 0`. The `False()` setting is
therefore unusable here — `(a × b) ⋅ b` is not bit-zero, so a point sitting
*exactly* on a vertex tests as off the arc, and every boundary case is missed.
And `True()` is unaffordable: it lifts to `Rational{BigInt}` whenever two arcs
share a great circle, ~11 KB per call, and two cells sharing an edge is the
common case in a tiling rather than a corner case.

A banded predicate is the right tier for a lightweight yes/no answer: it costs
nothing, and it decides the shared-vertex and shared-edge cases that a tiling
actually produces. Callers who need exact spherical topology want RelateNG.

## Ring semantics

A ring is read as the region it encloses, by even-odd crossing parity against a
definitionally exterior anchor ([`spherical_ring_encloses`](@ref)). That is the
spherical counterpart of the planar family's even-odd ray cast, and it is
winding-independent, so a CW and a CCW spelling of the same ring agree.
=#

"""
    SphericalRingPoints(ring)

A ring's vertices as `UnitSphericalPoint`s, converted on indexing rather than
up front, with any repeated closing vertex excluded from `length`.

Exists so that the spherical predicates can hand a ring to the
`UnitSpherical` ring primitives — which index into a vector of unit points —
without allocating that vector. Converting lazily costs one
`UnitSphereFromGeographic` per access; a ring walked several times in one
predicate call pays that more than once, which is the trade a prepared target
would remove.
"""
struct SphericalRingPoints{G} <: AbstractVector{UnitSphericalPoint{Float64}}
    ring::G
    n::Int
end

function SphericalRingPoints(ring)
    n = GI.npoint(ring)
    #= A closed ring repeats its first vertex; the primitives here take the
    closing edge as implied, so counting it would place a zero-length edge in
    the walk. =#
    if n > 1 && equals(GI.getpoint(ring, 1), GI.getpoint(ring, n))
        n -= 1
    end
    return SphericalRingPoints(ring, n)
end

Base.size(v::SphericalRingPoints) = (v.n,)
Base.IndexStyle(::Type{<:SphericalRingPoints}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::SphericalRingPoints, i::Int) =
    _spherical_kernel_point(GI.getpoint(v.ring, i))

#=
Whether `p` lies on the closed minor arc `a → b`.

Not `UnitSpherical.point_on_spherical_arc`, whose span test compares cosines
(`a⋅p ≥ a⋅b` and `b⋅p ≥ a⋅b`). That criterion is only sufficient for arcs
shorter than a quarter turn: for a long arc it also admits points on the far
side of `a`. The test arcs here run from a query point to a definitionally
exterior anchor and are therefore *close to a half turn* by construction, so the
failure is not exotic — it silently inverts crossing parity, and with it every
containment answer. (Measured: a ring vertex 9.9° west of the query counted as
lying on an arc running 175° east.)

The determinant span test `_on_arc_span_authority` is exact in form and correct
for any arc length; in `Float64` it is allocation-free. The great-circle gate in
front of it is the banded `spherical_orient`, keeping this in the same tolerance
regime as everything else in this file.
=#
@inline function _sph_on_arc(p, a, b)
    UnitSpherical.spherical_orient(a, b, p) == 0 || return false
    return _on_arc_span_authority(False(), p, a, b)
end

#=
Move the anchor off `q`'s antipode, if that is where it landed.

The parity walk runs a test arc from `q` to the anchor, which is undefined when
the two are antipodal — and they are *exactly* antipodal for the most natural
query there is: a ring's own centre, since the default anchor is the antipode of
the vertex mass. Without this, locating a cell's centroid inside its own cell
fails.

The anchor is a free choice of any point in the exterior region, so nudging it a
milliradian off the antipode is legitimate: a ring whose enclosed region falls
short of a hemisphere by more than that still has the nudged point outside. Rings
too close to a hemisphere for that margin are exactly the ones
`spherical_exterior_anchor` already declines to anchor at all.
=#
@inline function _nudge_anchor(z, q)
    z === nothing && return z
    dot(q, z) >= -1 + 1e-9 && return z
    # some direction not parallel to z, made orthogonal to it
    ref = abs(z[1]) < 0.9 ? UnitSphericalPoint(1.0, 0.0, 0.0) :
                            UnitSphericalPoint(0.0, 0.0, 1.0)
    t = ref - dot(ref, z) .* z
    nt = norm(t)
    nt == 0 && return z
    return UnitSphericalPoint(normalize(z + (1e-3 / nt) .* t))
end

#=
Locate `q` against the ring `v` (already unit points), returning `in`, `on` or
`out`.

Mirrors `rk_point_in_ring`, minus the prepared `SphericalKernelRing`: boundary
membership is settled first by a per-edge scan, then interior membership by
even-odd parity with boundary points already excluded.

`anchor` is passed in because a caller testing several points against one ring
(the hinge walk) can compute it once.
=#
function _usp_ring_orientation(
    m::Spherical, v, anchor, q;
    in::T = point_in, on::T = point_on, out::T = point_out, exact,
) where {T}
    n = length(v)
    n == 0 && return out
    @inbounds for j in 1:n
        _sph_on_arc(q, v[j], v[mod1(j + 1, n)]) && return on
    end
    # Fewer than three distinct vertices bound no area.
    n < 3 && return out
    enc = UnitSpherical.spherical_ring_encloses(v, n, q;
        anchor = _nudge_anchor(anchor, q),
        on_arc = Returns(false),      # boundary settled above
        on_test_arc = _sph_on_arc,    # the near-half-turn arc the default breaks on
    )
    enc === nothing && _throw_degenerate_ring_orientation(q)
    return enc ? in : out
end

@noinline _throw_degenerate_ring_orientation(q) = throw(ArgumentError(
    "the lightweight spherical predicates cannot locate the point $(q) against " *
    "this ring: its vertex mass is degenerate (a near-hemisphere or " *
    "vertex-symmetric ring), so no definitionally exterior anchor exists. Use " *
    "`relate_predicate(RelateNG(Spherical()), pred, a, b)`, which falls back to " *
    "a winding-consistent wedge bootstrap."))

function _point_filled_curve_orientation(
    m::Spherical, point, curve;
    in::T = point_in, on::T = point_on, out::T = point_out, exact,
) where {T}
    v = SphericalRingPoints(curve)
    q = _spherical_kernel_point(point)
    anchor = length(v) >= 3 ? UnitSpherical.spherical_exterior_anchor(v, length(v)) : nothing
    return _usp_ring_orientation(m, v, anchor, q; in, on, out, exact)
end

#=
`on` when `point` is one of the segment's endpoints, `in` when it lies strictly
between them along the great-circle arc, `out` otherwise — matching the planar
method's split.

A zero-length segment (`start == stop`, which real rings do carry near polar
corners) cannot contain a point other than that vertex, and reaching
`_sph_on_arc` with one is safe: it is a sign and determinant test throughout,
with no division by the segment direction.
=#
function _point_segment_orientation(
    m::Spherical, point, start, stop;
    in::T = point_in, on::T = point_on, out::T = point_out,
) where {T}
    q = _spherical_kernel_point(point)
    a = _spherical_kernel_point(start)
    b = _spherical_kernel_point(stop)
    (q == a || q == b) && return on
    return _sph_on_arc(q, a, b) ? in : out
end

#=
Map the symbolic arc classification onto the `(orientation, α, β)` the
processors read.

`α` and `β` locate the meeting along `(a1, a2)` and `(b1, b2)`, but callers only
ever ask whether either is `0` or `1` — whether the meeting is at a segment
endpoint. The endpoint incidences answer exactly that, so interior meetings
report `0.5`: a value that is neither endpoint, and is never compared for
anything else.
=#
function _seg_seg_orientation(m::Spherical, a1, a2, b1, b2; exact)
    ka1 = _spherical_kernel_point(a1)
    ka2 = _spherical_kernel_point(a2)
    kb1 = _spherical_kernel_point(b1)
    kb2 = _spherical_kernel_point(b2)
    orient, a0_on_b, a1_on_b, b0_on_a, b1_on_a =
        _sph_arc_arc_class(ka1, ka2, kb1, kb2)
    orient === line_out && return line_out, 0.5, 0.5
    α = a0_on_b ? 0.0 : (a1_on_b ? 1.0 : 0.5)
    β = b0_on_a ? 0.0 : (b1_on_a ? 1.0 : 0.5)
    return orient, α, β
end

#=
Classify two great-circle arcs, returning the `LineOrientation` and which of the
four endpoints lies on the other arc.

The shape of the test is the reduction `rk_classify_intersection` documents, run
on the banded predicates instead of the exact ones: four orientations place each
arc's endpoints against the other's great circle, endpoint incidences come from
arc membership, and a proper crossing is the strict straddle pattern in both
directions. No intersection coordinate is constructed.

Degenerate (zero-length) arcs are settled first. Real rings carry repeated
vertices — HEALPix rings do near polar corners — and every test below is a sign
or a dot product, so a zero-length arc produces no division and no NaN; it is
handled up front only because "the arc's normal" is meaningless for one, not
because it would misbehave.
=#
function _sph_arc_arc_class(a0, a1, b0, b1)
    adeg = a0 == a1
    bdeg = b0 == b1
    if adeg && bdeg
        same = a0 == b0
        return (same ? line_hinge : line_out), same, same, same, same
    elseif adeg
        on = _sph_on_arc(a0, b0, b1)
        return (on ? line_hinge : line_out), on, on, on && b0 == a0, on && b1 == a0
    elseif bdeg
        on = _sph_on_arc(b0, a0, a1)
        return (on ? line_hinge : line_out), on && a0 == b0, on && a1 == b0, on, on
    end

    sab0 = UnitSpherical.spherical_orient(a0, a1, b0)
    sab1 = UnitSpherical.spherical_orient(a0, a1, b1)
    sba0 = UnitSpherical.spherical_orient(b0, b1, a0)
    sba1 = UnitSpherical.spherical_orient(b0, b1, a1)

    a0_on_b = sba0 == 0 && _sph_on_arc(a0, b0, b1)
    a1_on_b = sba1 == 0 && _sph_on_arc(a1, b0, b1)
    b0_on_a = sab0 == 0 && _sph_on_arc(b0, a0, a1)
    b1_on_a = sab1 == 0 && _sph_on_arc(b1, a0, a1)

    if sab0 == 0 && sab1 == 0 && sba0 == 0 && sba1 == 0
        #= Same great circle. The arcs overlap iff some endpoint of one lies on
        the other; the overlap is a single point exactly when every incident
        endpoint is the same point, which is a hinge rather than an overlap. =#
        (a0_on_b || a1_on_b || b0_on_a || b1_on_a) ||
            return line_out, false, false, false, false
        shared = _sole_shared_point(a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
        return (shared ? line_hinge : line_over), a0_on_b, a1_on_b, b0_on_a, b1_on_a
    end

    (a0_on_b || a1_on_b || b0_on_a || b1_on_a) &&
        return line_hinge, a0_on_b, a1_on_b, b0_on_a, b1_on_a

    #= A transversal crossing. Each great circle separating the other's endpoints is
    necessary but *not* sufficient, because two great circles meet at an antipodal pair and
    a straddle in both directions does not say the two arcs reach the *same* one of them.

    Two small arcs on opposite sides of the sphere show it: `(-0.5,-0.5)→(0.5,-0.5)` and
    `(180.499,0.5)→(179.499,0.5)` are near-antipodal images of one another, so each lies
    almost on the other's great circle and straddles it, yet one arc contains one meeting
    point and the other arc contains its antipode. Under the two-way test they cross, and
    two disjoint cells on opposite sides of the globe report an intersection the size of a
    whole cell.

    The sign pattern below is S2's `SimpleCrossing`: with `acb = -sab0`, `bda = sab1`,
    `cbd = -sba1` and `dac = sba0`, it asks that `acb` agree in sign with all three, which
    pins both arcs to the same meeting point. =#
    if sab0 != 0 && sab1 == -sab0 && sba0 == -sab0 && sba1 == sab0
        return line_cross, false, false, false, false
    end
    return line_out, false, false, false, false
end

#= Whether every endpoint incidence names one and the same point — a collinear
abutment rather than an overlap of positive length. =#
@inline function _sole_shared_point(a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    p = a0_on_b ? a0 : (a1_on_b ? a1 : (b0_on_a ? b0 : b1))
    (!a0_on_b || a0 == p) && (!a1_on_b || a1 == p) &&
        (!b0_on_a || b0 == p) && (!b1_on_a || b1 == p)
end

#=
The unit point where arcs `(a0, a1)` and `(b0, b1)` cross, for a pair already
already classified as a proper crossing.

Two great circles meet at an antipodal pair `±x`; the classification guarantees
the crossing lies strictly inside both minor arcs, so exactly one of `±x` does,
and it is the one in the same hemisphere as either arc's midpoint direction.
That check is what makes choosing between the two candidates safe here — it is
sound *because* the crossing is already known to be proper, and would not be on
its own.

Returns `nothing` if the normals are parallel, which a proper crossing cannot
produce, but which guards the normalization regardless.

The two normals come from `robust_cross_product`, not from a plain `cross`. For
two vertices a cell edge apart the plain product is a difference of `O(1)` terms
whose result is `O(4e-6)` at HEALPix level 18, so it keeps only the bits the
cancellation leaves, and crossing two such normals compounds it. Measured
against a `BigFloat` reference on level-18 cells at longitude 120, the plain
form places the crossing 4.8e-7 of an edge length away (worst 2.8e-6); through
`robust_cross_product` — the same stable `(a-b) × (a+b)` form `spherical_orient`
uses, with the exact fallback behind it — that becomes 3e-11 (worst 9e-11), and
it still allocates nothing. Degree-scale inputs barely notice the difference,
which is exactly why this has to be measured at cell scale.
=#
@inline function _arc_crossing_point(a0, a1, b0, b1)
    x = cross(robust_cross_product(a0, a1), robust_cross_product(b0, b1))
    nx = norm(x)
    nx == 0 && return nothing
    u = x ./ nx
    return UnitSphericalPoint(dot(u, a0 + a1) < 0 ? -u : u)
end

#=
Walk the segment `l_start → l_end` across `curve`, splitting it at every point
where the two meet, and report whether the pieces run inside and/or outside the
filled curve.

The planar method materializes every intersection point, sorts them by distance
from the segment start, and tests each piece's midpoint. This does the same walk
without the vector or the sort: it repeatedly scans the ring for the split point
nearest the current position and steps to it.

`dot(A, ·)` decreases monotonically along a minor arc from `A`, so it orders
split points along the segment without any angle being computed; each step takes
the largest such value strictly below the current one. Strictly below is what
terminates the walk — coincident split points are visited once, and the position
advances every iteration.
=#
function _split_segment_interactions(
    m::Spherical, l_start, l_end, curve, in_curve, out_curve; exact,
)
    A = _spherical_kernel_point(l_start)
    B = _spherical_kernel_point(l_end)
    v = SphericalRingPoints(curve)
    n = length(v)
    (n == 0 || A == B) && return in_curve, out_curve
    anchor = n >= 3 ? UnitSpherical.spherical_exterior_anchor(v, n) : nothing

    t_end = dot(A, B)
    p_start = A
    #= The walk's position is tracked by its own ordinate, not by an idealized
    `1`. `dot(A, A)` is a rounded sum of three squares and lands an ulp or so
    below one, so seeding `t_start` with `one(t_end)` would leave the start
    point itself strictly ahead of the position: a ring vertex coincident with
    `A` — every shared edge has one — would then be taken as a split point and
    open a zero-length piece. =#
    t_start = dot(A, A)
    while true
        # the split point nearest the current position, if any is left
        best_t = t_end
        best_p = B
        found = false
        @inbounds for j in 1:n
            c0 = v[j]
            c1 = v[mod1(j + 1, n)]
            orient, _, _, b0_on_a, b1_on_a = _sph_arc_arc_class(A, B, c0, c1)
            orient === line_out && continue
            if orient === line_cross
                x = _arc_crossing_point(A, B, c0, c1)
                if x !== nothing
                    t = dot(A, x)
                    if t < t_start && t > best_t
                        best_t, best_p, found = t, x, true
                    end
                end
            else
                #= Any curve vertex lying on this segment splits it; the
                classifier's incidence flags name them, so no separate on-arc
                test (and no second tolerance) is needed. =#
                if b0_on_a
                    t = dot(A, c0)
                    if t < t_start && t > best_t
                        best_t, best_p, found = t, c0, true
                    end
                end
                if b1_on_a
                    t = dot(A, c1)
                    if t < t_start && t > best_t
                        best_t, best_p, found = t, c1, true
                    end
                end
            end
        end
        p_end = found ? best_p : B
        #= A piece with no extent says nothing about which side the segment
        runs, and cannot be asked: normalizing `p + p` moves the point by an
        ulp, and an ulp past a shared vertex is off the end of both arcs that
        meet there, so it would classify as outside the ring it is a vertex
        of. =#
        mid = p_start + p_end
        nm = norm(mid)
        if p_end != p_start && nm > 0
            mid_val = _usp_ring_orientation(m, v, anchor, UnitSphericalPoint(mid ./ nm); exact)
            if mid_val == point_in
                in_curve = true
            elseif mid_val == point_out
                out_curve = true
            end
        end
        found || break
        p_start = p_end
        t_start = best_t
    end
    return in_curve, out_curve
end

#=
A bounding cap for any geometry: centred on its normalized vertex mass, with the
radius that reaches its furthest vertex.

Sound as a *containing* bound because a spherical cap of radius under a quarter
turn is convex, so the minor arc between two vertices inside the cap stays
inside it. That is what makes this safe where a lon/lat box is not: a
great-circle arc bulges poleward of both its endpoints and can leave the box
that contains them, but it cannot leave the cap that contains them.

Two O(n) passes, no allocation. Returns `nothing` — no shortcut available —
when the vertex mass is degenerate or the cap would exceed a quarter turn, since
the convexity argument is what the soundness rests on.

The radius is built from the chord rather than `acos` of the dot product, which
loses all its precision for the small caps this is most useful on, and is then
padded: the shortcut may only ever be *more* reluctant to reject.
=#
_spherical_bounding_cap(geom) = _spherical_bounding_cap(GI.trait(geom), geom)

# `GI.getpoint` has no `PointTrait` method; a point is its own zero-radius cap.
_spherical_bounding_cap(::GI.PointTrait, geom) =
    UnitSpherical.SphericalCap(_spherical_kernel_point(geom), 1e-7)

function _spherical_bounding_cap(::GI.AbstractGeometryTrait, geom)
    sx = sy = sz = 0.0
    n = 0
    for p in GI.getpoint(geom)
        u = _spherical_kernel_point(p)
        sx += u[1]; sy += u[2]; sz += u[3]
        n += 1
    end
    n == 0 && return nothing
    nrm = sqrt(sx * sx + sy * sy + sz * sz)
    nrm < 1e-9 && return nothing  # vertices spread over a great circle
    c = UnitSphericalPoint(sx / nrm, sy / nrm, sz / nrm)
    maxchord2 = 0.0
    for p in GI.getpoint(geom)
        u = _spherical_kernel_point(p)
        dx = c[1] - u[1]; dy = c[2] - u[2]; dz = c[3] - u[3]
        ch2 = dx * dx + dy * dy + dz * dz
        ch2 > maxchord2 && (maxchord2 = ch2)
    end
    maxchord2 >= 2.0 && return nothing  # radius past a quarter turn
    r = 2 * asin(sqrt(maxchord2) / 2)
    return UnitSpherical.SphericalCap(c, r + 1e-7 + 8 * eps(r))
end

#=
The planar extent shortcut is unsound on the sphere, so the spherical one
rejects on bounding caps instead.

Two lon/lat boxes being disjoint does not make the geometries disjoint — a
great-circle edge can reach outside the latitude range of its own endpoints — so
inheriting the planar test would answer "disjoint" for geometries that meet.
Bounding caps carry the same cheap-rejection benefit with a bound that holds on
the sphere. The disposition of the answer, once the two are known disjoint, is
the planar method's unchanged.

This matters out of proportion to its size for the workload this path exists
for: testing many thousands of cells against one target, where nearly every
candidate is far away and never needed the full walk.
=#
@inline function _maybe_skip_disjoint_extents(::Spherical, a, b;
    in_allow, on_allow, out_allow,
    in_require, on_require, out_require,
    kw...
)
    capa = _spherical_bounding_cap(a)
    capa === nothing && return (false, false)
    capb = _spherical_bounding_cap(b)
    capb === nothing && return (false, false)
    UnitSpherical._disjoint(capa, capb) || return (false, false)
    return if out_allow
        (in_require || on_require) ? (true, false) : (true, true)
    else
        # points not allowed in the exterior, but the geometries are disjoint
        (true, false)
    end
end
