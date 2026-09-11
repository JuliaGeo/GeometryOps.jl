# # Precompile workload
#
#=
Precompile RelateNG's shared kernel engine and native geometry ingestion paths. Predicate
types specialize the topology computer, so each predicate needs a call.

OverlayNG uses one operation value to compile the shared driver. Each manifold needs separate
instances; geometry shapes cover ingestion, line building, and mixed-point paths.

Cache only the default output point type. Spherical Foster-Hormann workloads cover lon/lat and
xyz input, polygon tracing, and area measurement.
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

    #-- Spherical Foster–Hormann has separate lon/lat ingestion and Cartesian
    #-- point paths. Use overlapping polygons so both trace proper crossings.
    _pc_fh_ll = map((_pc_poly1, _pc_poly2)) do poly
        GI.Polygon([collect(GI.getpoint(GI.getexterior(poly)))])
    end
    _pc_fh_xyz = map(_pc_fh_ll) do poly
        apply(UnitSpherical.UnitSphericalPoint, GI.PointTrait(), poly)
    end

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

        #-- Compile area, line, and mixed-point paths on each manifold.
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

        #-- Compile Float64 polygon tracing and area measurement for both input representations.
        fh = FosterHormannClipping(Spherical())
        for (a, b) in (_pc_fh_ll, _pc_fh_xyz)
            cache = FosterHormannCache(fh)
            intersection_area(fh, a, b; cache)
            intersection_area(fh, a, b)
            intersection(fh, a, b; target = GI.PolygonTrait())
        end
    end
end
