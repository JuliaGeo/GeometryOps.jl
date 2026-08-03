# # Precompile workload
#
#=
First-call latency of the RelateNG predicates is dominated by inferring the
engine (topology computer, edge intersector, tree traversals) plus one thin
per-geometry-type outer layer (`RelateGeometry` construction, extraction).
The engine core is typed on kernel-level types only — see the opaque
geometry references in `TopologyComputer` / `RelateSegmentString` — so one
workload run caches it for *every* input geometry type; the outer layer is
exercised here for the native geometry types (`GO.tuples` output wrapped in
`GI.Wrappers`, which is also what the tests feed). This matters most on
Julia 1.12, where inference of these instances is several times slower than
on 1.11.

The OverlayNG block below is the same shape of workload for the overlay engine,
and is sized by measurement, not by guesswork — `benchmarks/overlayng_ttfx.jl`
is the probe, and its recorded before/after tables are the justification for
every line of it. Three facts from that probe determine the shape:

1. The op is a *value*, not a type (`_overlay_ng(m, op::_OverlayOpCode, a, b)`),
   so one call caches the driver for all four operations: with only
   `intersection` in the workload, first-call `union`/`difference`/
   `symdifference` on the same input types drops from ~1.1 s to ~2 ms (the
   residual is the thin op wrapper). Precompiling the other three would buy
   ~6 ms and is not worth an instance.
2. The arrangement, graph and builders are typed on the kernel point type only
   (`NodedArrangement{P}`, `OverlayGraph`), exactly like RelateNG's engine, so
   one call per manifold caches the whole engine core for every input geometry
   type. What remains per input type is the ingest layer, which is why the
   `mpoly` and `line` shapes below are nearly free (~0.4 MB for the two of them
   on the plane, against 2.0 MB for the first `poly x poly` call).
3. Manifolds do *not* share instances, and the point path
   (`_overlay_mixed_points`) never reaches the arrangement, so those are the
   two genuinely separate cost centres — and both pay for themselves; see the
   ledger in `benchmarks/overlayng_ttfx.jl`.

Cost measured on this machine (Apple M4 Pro, Julia 1.12.6): the OverlayNG block
adds ~4.4 s to the package precompile (14.6 s → 19.0 s) and ~8.9 MB to the
pkgimage (24.0 → 32.8 MB, which shows up as ~40 ms of extra package load), and
takes the summed first call over a 32-instance op × shape × manifold matrix from
45.0 s to 0.44 s.
=#

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    _pc_ring(pts) = GI.LinearRing(pts)
    _pc_poly1 = GI.Polygon([_pc_ring([(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0), (0.0, 0.0)])])
    _pc_poly2 = GI.Polygon([_pc_ring([(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0), (2.0, 2.0)])])
    _pc_mpoly = GI.MultiPolygon([_pc_poly1, _pc_poly2])
    _pc_line = GI.LineString([(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])
    _pc_mline = GI.MultiLineString([_pc_line, GI.LineString([(0.0, 1.0), (2.0, 1.0)])])
    _pc_pt = GI.Point((1.0, 1.0))
    _pc_geoms = (_pc_poly1, _pc_mpoly, _pc_line, _pc_mline, _pc_pt)
    #-- overlay contracts on valid input, so its multipolygon needs disjoint
    #-- components (`_pc_mpoly` above deliberately overlaps, for `relate`)
    _pc_poly3 = GI.Polygon([_pc_ring([(6.0, 0.0), (8.0, 0.0), (8.0, 2.0), (6.0, 2.0), (6.0, 0.0)])])
    _pc_mpoly_d = GI.MultiPolygon([_pc_poly1, _pc_poly3])
    _pc_ovl_line = GI.LineString([(-1.0, 1.0), (1.5, 1.5), (4.0, 1.0)])

    @compile_workload begin
        alg = RelateNG()
        #-- every predicate re-specializes the topology computer on its
        #-- predicate type; one polygon-pair call each caches the engine
        for f in (intersects, disjoint, contains, within, covers,
                coveredby, crosses, overlaps, touches, equals)
            f(alg, _pc_poly1, _pc_poly2)
        end
        #-- the per-geometry-type outer layer (RelateGeometry construction,
        #-- extraction, point location), over the native type combinations
        for a in _pc_geoms, b in _pc_geoms
            relate(alg, a, b)
        end
        #-- prepared mode
        prep = prepare(alg, _pc_poly1)
        relate(prep, _pc_poly2)
        relate(prep, _pc_pt)

        #-- OverlayNG: one op (the op code is a value — see the note above)
        #-- over the four input shapes that reach different code, on each
        #-- manifold. Area x area runs the whole arrangement; line x area adds
        #-- the line builder; point x area takes the separate mixed-points path.
        for m in (Planar(), Spherical())
            ovl = OverlayNG(m)
            intersection(ovl, _pc_poly1, _pc_poly2)
            intersection(ovl, _pc_mpoly_d, _pc_poly2)
            intersection(ovl, _pc_ovl_line, _pc_poly1)
            intersection(ovl, _pc_pt, _pc_poly1)
            #-- a target is a singleton type, so it specializes the driver and
            #-- the extractor afresh; the areal one is the case worth caching
            intersection(ovl, _pc_poly1, _pc_poly2; target = GI.MultiPolygonTrait())
        end
    end
end
