# # Antimeridian splitting
export antimeridian_split

#=
## What is the antimeridian problem?

A longitude/latitude polygon that spans the ±180° meridian (the *antimeridian*,
or *seam*) cannot be drawn as a single planar ring: a vertex at 179°E and the
next at 179°W are neighbours on the globe, but 358° apart on the page. Renderers,
GeoJSON (RFC 7946 §3.1.9), and most planar clippers all assume consecutive
vertices are close in longitude, so a seam-crossing polygon paints a spurious
band right across the map.

The fix is to *cut* the polygon along the seam into pieces that each live within
one 360°-wide longitude branch. This is an **encoding repair, not a geometry
correction**: the spherical region is unchanged — we only re-express it as a
`MultiPolygon` whose parts avoid the seam. Because the trait of the output
differs from the input (a `Polygon` can become a two-part `MultiPolygon`), this
is a *transformation* that returns a new geometry, not a
`GeometryCorrection` that repairs in place.

## Why the pole row?

A piece that *encloses a pole* is bounded, on the seam side, by an edge that runs
**along the pole** — the polygon walks up the seam to (say) +180°, crosses the
pole, and comes back down at −180°. On the sphere that pole edge is a single
point (a zero-length geodesic: both corners are the same kernel point), so it
carries no area. But pseudocylindrical projections (plate carrée, Robinson, …)
map the pole to a *line*, and a two-vertex pole edge would draw as a single
straight segment cutting across the top/bottom of the map. Resampling the pole
edge into a monotone row of constant-latitude points (`pole_spacing`) makes that
line trace the projected pole correctly. The samples are collinear no-ops in
plate carrée and have exactly zero spherical area, so they never change the
region — they exist only so the projected outline is right.

## Composition with `segmentize`

Run [`antimeridian_split`](@ref) **after** [`segmentize`](@ref), never before.
Spherical `segmentize` of an on-meridian edge emits interpolated points whose
`atan2` longitude is +180 regardless of which branch the edge belongs to
(signed-zero normalisation), which would corrupt a −180-branch piece. Densify
first, split second.

## Arbitrary seam and rotated pole

`antimeridian = λ` cuts along the meridian `λ` instead of ±180. The input
coordinates are left **untouched** — only the cutting curve and the output
longitude encoding move — so exactness is fully preserved. `north_pole =
(λp, φp)` cuts along a seam through a *rotated* pole: the input is rotated into
the frame whose north pole sits at geographic `(λp, φp)` (the CF
`rotated_latitude_longitude` / PROJ `+proj=ob_tran` convention, with
`north_pole_grid_longitude = 0`), the split runs in that frame, and **the output
coordinates are in the rotated frame** — that is the point of the feature.
Mapping back would re-merge the seam and produce degenerate slit lips. The
rotation is the *only* step that moves a coordinate; everything else is the exact
arrangement operating on the (rotated) input.

## Example

```@example antimeridian
import GeometryOps as GO
import GeoInterface as GI

# a box straddling the ±180° seam (170°E … 170°W)
box = GI.Polygon([[(170.0, 40.0), (-170.0, 40.0), (-170.0, 50.0), (170.0, 50.0), (170.0, 40.0)]])
mp = GO.antimeridian_split(box)
GI.ngeom(mp)   # 2 pieces, one per branch
```
=#

# The split is inherently spherical; these are the manifold/exactness context the
# arrangement runs under.
const _AM_MANIFOLD = Spherical()
const _AM_EXACT = True()

# lat within this of ±90 is treated as the pole; lon within this of the seam
# meridian (mod 360) is treated as on-seam.
const _AM_POLE_TOL = 1e-9
const _AM_SEAM_TOL = 1e-9

"""
    antimeridian_split(geom; antimeridian = 180.0, north_pole = nothing, pole_spacing = 5.0)

Split a lon/lat `Polygon` or `MultiPolygon` at the antimeridian, returning a
`GI.MultiPolygon` whose pieces each stay within one 360°-wide longitude branch
(none crosses the seam). The spherical region is preserved exactly — this is an
encoding repair, not a geometry change.

## Keyword arguments

- `antimeridian = 180.0`: the seam longitude `λ`. The cut runs along the meridian
  `λ` (normalised to `λn ∈ (-180, 180]`); the default reproduces the ±180°
  antimeridian. Emitted longitudes lie in the closed branch `[λn - 360, λn]`: the
  two seam lips are exactly `λn` (west-side pieces) and `λn - 360` (east-side
  pieces), and every non-seam vertex is strictly interior, in `(λn - 360, λn)`.
  At the default seam this is the usual `[-180, 180]` with lips at +180 / −180.
- `north_pole = nothing`: when set to `(λp, φp)`, cut along a seam through the
  rotated pole at geographic `(λp, φp)` (CF `rotated_latitude_longitude` /
  `+proj=ob_tran`, `north_pole_grid_longitude = 0`). **The returned coordinates
  are in the rotated frame**, not geographic — this is deliberate.
- `pole_spacing = 5.0`: maximum longitude step (degrees) of the constant-latitude
  row emitted along a pole edge of a pole-enclosing piece. `nothing` emits only
  the two branch corners. The corners are never optional — they are the
  topological product of the face walk; only the infill between them is
  controlled here.

Only `PolygonTrait` and `MultiPolygonTrait` inputs are supported; other traits
throw an `ArgumentError` (LineString support is future work).
"""
function antimeridian_split(geom; antimeridian = 180.0, north_pole = nothing, pole_spacing = 5.0)
    t = GI.trait(geom)
    (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) || throw(ArgumentError(
        "antimeridian_split supports PolygonTrait and MultiPolygonTrait inputs; " *
        "got $(t === nothing ? typeof(geom) : typeof(t)). (LineString support is future work.)"))
    λn = _normalize_seam(Float64(antimeridian))
    work = north_pole === nothing ? geom :
           _rotate_to_pole(geom, Float64(north_pole[1]), Float64(north_pole[2]))
    pieces = Vector{Vector{Vector{Tuple{Float64, Float64}}}}()
    _each_polygon(work) do poly
        _split_polygon!(pieces, poly, λn; pole_spacing)
    end
    return GI.MultiPolygon(pieces)
end

# Normalise a seam longitude into `(-180, 180]` (so `±180` collapse to `180`).
function _normalize_seam(λ::Float64)
    m = rem(λ, 360.0, RoundNearest)      # [-180, 180]
    return m == -180.0 ? 180.0 : m
end

# The pole-to-pole cutting arc at meridian `λn`, as two quarter-great-circle
# segments (a single antipodal segment has no unique arc and is rejected at
# ingest).
_meridian_arc(λn) = GI.LineString([(λn, -90.0), (λn, 0.0), (λn, 90.0)])

# Iterate the polygon parts of a Polygon/MultiPolygon.
function _each_polygon(f::F, geom) where {F}
    if GI.trait(geom) isa GI.PolygonTrait
        f(geom)
    else
        for poly in GI.getgeom(geom)
            f(poly)
        end
    end
    return nothing
end

# The exterior ring followed by the hole rings of a polygon.
_rings_of(poly) = (GI.getexterior(poly), GI.gethole(poly)...)

_ring_coords(ring) = Tuple{Float64, Float64}[(Float64(GI.x(p)), Float64(GI.y(p)))
                                             for p in GI.getpoint(ring)]

# ## Rotated-pole frame (tier 2)
#
# The only coordinate-moving step in the whole transformation. Rotate a
# geographic unit vector into the frame whose north pole is at geographic
# `(λp, φp)` (CF `rotated_latitude_longitude`, `north_pole_grid_longitude = 0`):
# `R = R_y(φp − 90°) · R_z(−λp)`, so `R · P = (0, 0, 1)` for `P` the unit vector
# of `(λp, φp)`. `R` is a proper rotation (det = +1), so ring orientation and
# area are preserved.
function _pole_rotation_matrix(λp::Float64, φp::Float64)
    sλ, cλ = sincosd(λp)
    sφ, cφ = sincosd(φp)
    return @SMatrix [ sφ*cλ   sφ*sλ  -cφ ;
                     -sλ      cλ      0.0 ;
                      cφ*cλ   cφ*sλ   sφ ]
end

function _rotate_to_pole(geom, λp::Float64, φp::Float64)
    R = _pole_rotation_matrix(λp, φp)
    fromgeo = UnitSpherical.UnitSphereFromGeographic()
    togeo = UnitSpherical.GeographicFromUnitSphere()
    return apply(GI.PointTrait(), geom) do p
        v = fromgeo((Float64(GI.x(p)), Float64(GI.y(p))))
        togeo(R * v)
    end
end

# ## Per-part split (tier 1, frame-agnostic)
#
# Each polygon part is split independently; non-interacting parts pass through
# untouched (fast path).
function _split_polygon!(pieces, poly, λn; pole_spacing)
    if !_may_cross_seam(poly, λn)
        push!(pieces, [_ring_coords(r) for r in _rings_of(poly)])
        return nothing
    end
    arc = _meridian_arc(λn)
    ssa = _overlay_segstrings(_AM_MANIFOLD, poly, true; exact = _AM_EXACT)
    ssb = _overlay_segstrings(_AM_MANIFOLD, arc, false; exact = _AM_EXACT)
    arr = NodedArrangement(_AM_MANIFOLD, ssa, ssb; exact = _AM_EXACT)
    g = OverlayGraph(_AM_MANIFOLD, arr; exact = _AM_EXACT)
    input = _OverlayInput(_AM_MANIFOLD, poly, arc, 2, 1, _AM_EXACT, false, false, nothing, nothing)
    _compute_labelling!(g, input)
    ctx = _build_faces(_AM_MANIFOLD, g; exact = _AM_EXACT)

    # keep the faces whose A-side (polygon) location is interior; CW rings are
    # shells, CCW rings are cavities assigned to their shells by the shared
    # containment machinery.
    for er in 1:length(ctx.edge_rings)
        _face_ring_location(ctx, er, 0) == LOC_INTERIOR || continue
        ring = ctx.edge_rings[er]
        ring.is_hole ? push!(ctx.free_hole_list, Int32(er)) :
                       push!(ctx.shell_list, Int32(er))
    end
    _place_free_holes!(ctx)

    for sh in ctx.shell_list
        rings = [_emit_ring(ctx, g, Int(sh), λn; pole_spacing)]
        for h in ctx.edge_rings[sh].holes
            push!(rings, _emit_ring(ctx, g, Int(h), λn; pole_spacing))
        end
        push!(pieces, rings)
    end
    return nothing
end

# ## Fast path — does a part interact with the seam meridian?
#
# A part interacts iff some ring edge's minor great-circle arc can reach the seam
# half-plane (meridian `λn`, i.e. `{v · w = 0, v · u > 0}` for the equatorial
# unit vectors `u` toward `λn` and `w` toward `λn + 90°`). Because an open
# hemisphere is geodesically convex, a minor arc whose endpoints are both
# strictly on one side of the seam great circle (`v · w`), or both strictly in
# the far hemisphere (`v · u < 0`, toward `λn − 180`), stays there and cannot
# touch the seam. Anything else routes through the full machinery — false
# positives only cost speed, never correctness. This is a spherical (kernel
# point) test, not a planar longitude heuristic: it detects pole-enclosing rings
# (which must wind 360° in longitude, hence cross the seam) and seam-crossing
# near-pole edges that a `|Δlon|` threshold misses.
function _may_cross_seam(poly, λn)
    su, cu = sincosd(λn)
    ux, uy = cu, su          # u = (cos λn, sin λn, 0), toward the seam meridian
    wx, wy = -su, cu         # w = (-sin λn, cos λn, 0), normal to the seam plane
    fromgeo = UnitSpherical.UnitSphereFromGeographic()
    for ring in _rings_of(poly)
        du_prev = dw_prev = 0.0
        have_prev = false
        first_du = first_dw = 0.0
        for p in GI.getpoint(ring)
            v = fromgeo((Float64(GI.x(p)), Float64(GI.y(p))))
            du = v[1] * ux + v[2] * uy
            dw = v[1] * wx + v[2] * wy
            if have_prev
                _arc_may_cross_seam(du_prev, dw_prev, du, dw) && return true
            else
                first_du, first_dw = du, dw
                have_prev = true
            end
            du_prev, dw_prev = du, dw
        end
        # closing edge (harmless when the ring is already closed)
        have_prev && _arc_may_cross_seam(du_prev, dw_prev, first_du, first_dw) && return true
    end
    return false
end

@inline function _arc_may_cross_seam(dua, dwa, dub, dwb)
    (dwa > 0 && dwb > 0) && return false   # arc stays on one side of the seam plane
    (dwa < 0 && dwb < 0) && return false
    (dua < 0 && dub < 0) && return false   # arc stays in the far (λn−180) hemisphere
    return true
end

# ## Emission — longitude-branch fixing and pole-pair insertion
#
# The face lies on the RIGHT of every ring half-edge, so a southward seam
# half-edge has its face to the west (branch `λn`) and a northward one to the
# east (branch `λn − 360`). Non-seam vertices are wrapped into the output branch
# `(λn − 360, λn]`; seam vertices take the branch of an incident cutting-arc
# edge (with a neighbour-side fallback for isolated touches); a pole vertex where
# the seam turns around emits as the two-branch pole pair, optionally resampled.
function _emit_ring(ctx, g, er::Integer, λn::Float64; pole_spacing = 5.0)
    edges = g.edges
    ring = ctx.edge_rings[er]
    es = Int32[]
    e = ring.start_edge
    while true
        push!(es, e)
        e = oe_next_result(edges, e)
        e == ring.start_edge && break
    end
    n = length(es)
    raw = [node_point(g.arr, he_origin(edges, es[j])) for j in 1:n]
    pts = Tuple{Float64, Float64}[]
    for j in 1:n
        (lon, lat) = raw[j]
        e_out = es[j]
        e_in = es[mod1(j - 1, n)]
        if _is_pole(lat)
            b_in = _seam_branch(g, e_in, λn)
            b_out = _seam_branch(g, e_out, λn)
            b1 = b_in === nothing ? b_out : b_in
            b2 = b_out === nothing ? b_in : b_out
            polelat = lat > 0 ? 90.0 : -90.0
            if b1 === nothing
                # pole vertex with no seam context (the ring touches the pole
                # without enclosing it): use the nearest non-pole neighbour's
                # branch so the emitted ring has no gratuitous longitude jump.
                push!(pts, (_neighbor_lon(raw, j, λn), polelat))
            elseif b1 == b2
                push!(pts, (b1, polelat))
            else
                push!(pts, (b1, polelat))                 # entry corner
                if pole_spacing !== nothing
                    # densified pole row: strictly monotone longitude from b1 to
                    # b2 through the branch interior (|b2 − b1| == 360), step ≤
                    # pole_spacing, both corners bit-exact.
                    ns = max(2, ceil(Int, 360.0 / Float64(pole_spacing)))
                    for k in 1:(ns - 1)
                        push!(pts, (b1 + (b2 - b1) * k / ns, polelat))
                    end
                end
                push!(pts, (b2, polelat))                 # exit corner
            end
        elseif _on_seam(lon, λn)
            b = _seam_branch(g, e_in, λn)
            b === nothing && (b = _seam_branch(g, e_out, λn))
            b === nothing && (b = _neighbor_lip(raw, j, λn))   # isolated touch
            push!(pts, (b, lat))
        else
            push!(pts, (_wrap_lon(lon, λn), lat))
        end
    end
    push!(pts, pts[1])
    return pts
end

# The seam branch of a half-edge that lies on the cutting arc (`nothing` if it is
# not part of the arc): south-running ⇒ west lip `λn`, north-running ⇒ east lip
# `λn − 360`.
function _seam_branch(g, e, λn::Float64)
    lbl = oe_label(g.edges, e)
    is_known(lbl, 1) || return nothing        # not part of the cutting arc (input B)
    (_, lat0) = node_point(g.arr, he_origin(g.edges, e))
    (_, lat1) = node_point(g.arr, he_dest(g.edges, e))
    return lat1 < lat0 ? λn : λn - 360.0
end

@inline _is_pole(lat) = abs(lat) >= 90.0 - _AM_POLE_TOL
@inline _on_seam(lon, λn) = abs(rem(lon - λn, 360.0, RoundNearest)) < _AM_SEAM_TOL

# Wrap a longitude into the output branch `(λn − 360, λn]`. Assumes `lon` is in
# `(-180, 180]` (as realised by `node_point`), so at the default seam `λn = 180`
# every in-range longitude is returned bit-unchanged.
@inline _wrap_lon(lon, λn) = lon > λn ? lon - 360.0 : lon

# The seam lip nearest a vertex's non-seam, non-pole neighbour (for an isolated
# seam touch): whichever of `λn` / `λn − 360` the wrapped neighbour longitude is
# closer to, split at the branch midpoint `λn − 180`.
function _neighbor_lip(raw, j, λn::Float64)
    n = length(raw)
    for k in 1:(n - 1)
        for idx in (mod1(j + k, n), mod1(j - k, n))
            lon = raw[idx][1]
            if !_on_seam(lon, λn) && !_is_pole(raw[idx][2])
                return _wrap_lon(lon, λn) > λn - 180.0 ? λn : λn - 360.0
            end
        end
    end
    return λn   # fully-degenerate ring; arbitrary
end

# The nearest non-pole ring vertex's wrapped longitude (for a pole vertex with no
# seam context), preferring the previous vertex so the incoming edge stays vertical.
function _neighbor_lon(raw, j, λn::Float64)
    n = length(raw)
    for k in 1:(n - 1)
        for idx in (mod1(j - k, n), mod1(j + k, n))
            _is_pole(raw[idx][2]) || return _wrap_lon(raw[idx][1], λn)
        end
    end
    return _wrap_lon(0.0, λn)
end
