# # Spherical leaves of the lightweight geometry-relation processors

#=
Spherical geometry methods for the manifold-independent DE-9IM allow/require processors in
`geom_geom_processors.jl`:

| leaf | question |
|:--|:--|
| `_point_segment_orientation` | is a point on a segment, and is it an endpoint? |
| `_point_filled_curve_orientation` | is a point in, on, or out of a filled ring? |
| `_seg_seg_orientation` | how do two segments meet? |
| `_split_segment_interactions` | does a hinging segment run inside, outside, or both? |

`SphericalRingPoints` converts vertices lazily. Segment classification uses endpoint
incidences without constructing intersections; segment splitting constructs only confirmed
proper crossings.

## Which tolerance regime, and why not the exact kernel's

The default uses `spherical_orient` with a `16 * eps` relative band (about 3.5e-15 radians)
and a determinant span test for all minor-arc lengths.

The band handles rounded shared-vertex and shared-edge contacts without exact arithmetic. Use
RelateNG for exact spherical topology.

## Ring semantics

Rings use winding-independent even-odd parity against an exterior anchor; see
[`spherical_ring_encloses`](@ref).
=#

"""
    SphericalRingPoints(ring)

Index a ring as `UnitSphericalPoint`s, excluding a repeated closing vertex from `length`.
Conversion occurs on every access without allocating a vertex vector.
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
Test membership on the closed minor arc `a → b`. A cosine span test also accepts points behind
an endpoint when the arc exceeds a quarter turn.

`_on_arc_span_authority` uses determinants valid for all minor arcs. The great-circle gate
uses the orientation selected by `exact`, defaulting to banded `spherical_orient`.
=#
@inline function _sph_on_arc(p, a, b, exact = False())
    _spherical_orient_for(booltype(exact))(a, b, p) == 0 || return false
    return _on_arc_span_authority(False(), p, a, b)
end

#=
Move an anchor away from `q`'s antipode, where the test arc is undefined. The default anchor
is antipodal to a query at the vertex-mass center.

A milliradian displacement remains exterior when the enclosed region is sufficiently smaller
than a hemisphere. `spherical_exterior_anchor` rejects near-degenerate vertex masses.
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
Locate `q` against unit-point ring `v`, returning `in`, `on`, or `out`. Test boundary
membership first, then even-odd interior parity.

Accept a precomputed `anchor` for repeated queries against the ring.
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
Return `on` for an endpoint, `in` for the arc interior, or `out` otherwise. A zero-length
segment contains only its endpoint; the determinant tests require no division.
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
Return `(orientation, α, β)` from symbolic arc classification. Endpoint incidences map to 0 or
1; interior meetings use 0.5 because callers test only endpoint equality.
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
Return the `LineOrientation` and four endpoint-incidence flags. Four orientations and arc
membership determine the result without constructing an intersection.

The default uses banded `spherical_orient`. Clipping requests exact signs because the band can
classify cell-scale crossings as hinges.

Handle zero-length arcs first because they have no great-circle normal.
=#
function _sph_arc_arc_class(a0, a1, b0, b1, exact = False())
    orient_of = _spherical_orient_for(booltype(exact))
    adeg = a0 == a1
    bdeg = b0 == b1
    if adeg && bdeg
        same = a0 == b0
        return (same ? line_hinge : line_out), same, same, same, same
    elseif adeg
        on = _sph_on_arc(a0, b0, b1, exact)
        return (on ? line_hinge : line_out), on, on, on && b0 == a0, on && b1 == a0
    elseif bdeg
        on = _sph_on_arc(b0, a0, a1, exact)
        return (on ? line_hinge : line_out), on && a0 == b0, on && a1 == b0, on, on
    end

    sab0 = orient_of(a0, a1, b0)
    sab1 = orient_of(a0, a1, b1)
    sba0 = orient_of(b0, b1, a0)
    sba1 = orient_of(b0, b1, a1)

    a0_on_b = sba0 == 0 && _sph_on_arc(a0, b0, b1, exact)
    a1_on_b = sba1 == 0 && _sph_on_arc(a1, b0, b1, exact)
    b0_on_a = sab0 == 0 && _sph_on_arc(b0, a0, a1, exact)
    b1_on_a = sab1 == 0 && _sph_on_arc(b1, a0, a1, exact)

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

    #= Both arcs must contain the same member of the antipodal intersection pair.
    S2's `SimpleCrossing` requires matching signs for `-sab0`, `sab1`, `-sba1`,
    and `sba0`; straddling each circle alone is insufficient. =#
    if sab0 != 0 && sab1 == -sab0 && sba0 == -sab0 && sba1 == sab0
        return line_cross, false, false, false, false
    end
    return line_out, false, false, false, false
end

# Test whether all endpoint incidences name one shared point.
@inline function _sole_shared_point(a0, a1, b0, b1, a0_on_b, a1_on_b, b0_on_a, b1_on_a)
    p = a0_on_b ? a0 : (a1_on_b ? a1 : (b0_on_a ? b0 : b1))
    (!a0_on_b || a0 == p) && (!a1_on_b || a1 == p) &&
        (!b0_on_a || b0 == p) && (!b1_on_a || b1 == p)
end

#=
Construct the unit intersection of arcs already classified as a proper crossing. Select the
antipodal candidate in the same hemisphere as either arc midpoint.

Return `nothing` for parallel normals. Use `robust_cross_product` to limit cancellation for
nearby endpoints.
=#
@inline function _arc_crossing_point(a0, a1, b0, b1)
    x = cross(robust_cross_product(a0, a1), robust_cross_product(b0, b1))
    nx = norm(x)
    nx == 0 && return nothing
    u = x ./ nx
    return UnitSphericalPoint(dot(u, a0 + a1) < 0 ? -u : u)
end

#=
Split `l_start → l_end` at contacts with `curve` and classify the pieces as inside or outside.
Repeated scans select the next split point without a vector or sort.

`dot(A, ·)` decreases along the minor arc. Select the largest value strictly below the current
one so coincident contacts are visited once and the walk advances.
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
    #= Start at `dot(A, A)`: rounding can place it below 1. Using 1 would revisit
    a coincident start vertex and create a zero-length piece. =#
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
        #= Skip zero-length pieces. Normalizing their midpoint can move a shared
        vertex off both incident arcs and misclassify it as exterior. =#
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
Bound the geometry by a cap centered on its normalized vertex mass and reaching its furthest
vertex. Caps smaller than a quarter turn contain the minor arcs between their vertices.

Return `nothing` for degenerate vertex mass or a cap reaching a quarter turn. Two passes take
O(n) time without allocation.

Compute the radius from chord length for small-angle precision, then pad it to keep rejection
conservative.
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
Reject disjoint spherical bounding caps. Endpoint lon/lat boxes are insufficient because
great-circle arcs can leave them. Once disjointness is established, use the planar flag rules.
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
