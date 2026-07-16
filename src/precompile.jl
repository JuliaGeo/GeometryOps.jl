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
    end
end
