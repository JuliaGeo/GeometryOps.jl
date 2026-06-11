# Tests for the RelateNG topology predicate framework
# (port of JTS TopologyPredicate / BasicPredicate / IMPredicate,
# plus the `intersects` and `disjoint` BasicPredicate kinds from
# RelatePredicate.java).

using Test
import GeometryOps as GO
import Extents

@testset "tri-state and intersects predicate" begin
    p = GO.pred_intersects()
    @test GO.predicate_name(p) == "intersects"
    @test !GO.is_known(p)
    @test GO.require_self_noding(typeof(p)) == false
    @test GO.require_interaction(typeof(p)) == true
    @test GO.require_covers(typeof(p), true) == false
    @test GO.require_exterior_check(typeof(p), true) == false
    @test GO.require_exterior_check(typeof(p), false) == false
    # interior/interior intersection determines intersects=true immediately
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # exterior-only updates never determine it; finish! defaults to false
    q = GO.pred_intersects()
    GO.update_dim!(q, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L)
    @test !GO.is_known(q)
    GO.finish!(q)
    @test GO.is_known(q) && GO.predicate_value(q) == false
    # disjoint envelopes determine intersects=false at init_bounds!
    r = GO.pred_intersects()
    GO.init_bounds!(r, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(5.0, 6.0), Y=(5.0, 6.0)))
    @test GO.is_known(r) && GO.predicate_value(r) == false
    # overlapping envelopes leave it unknown
    s = GO.pred_intersects()
    GO.init_bounds!(s, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(0.5, 6.0), Y=(0.5, 6.0)))
    @test !GO.is_known(s)
end

@testset "disjoint predicate" begin
    p = GO.pred_disjoint()
    @test GO.predicate_name(p) == "disjoint"
    @test GO.require_self_noding(typeof(p)) == false
    @test GO.require_interaction(typeof(p)) == false
    @test GO.require_exterior_check(typeof(p), true) == false
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    q = GO.pred_disjoint()
    GO.finish!(q)
    @test GO.is_known(q) && GO.predicate_value(q) == true
    # disjoint envelopes determine disjoint=true at init_bounds!
    r = GO.pred_disjoint()
    GO.init_bounds!(r, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(5.0, 6.0), Y=(5.0, 6.0)))
    @test GO.is_known(r) && GO.predicate_value(r) == true
end

@testset "tri-state value is sticky" begin
    # JTS BasicPredicate.setValue does not change an already-known value
    p = GO.pred_intersects()
    GO.update_dim!(p, GO.LOC_BOUNDARY, GO.LOC_BOUNDARY, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == true
    GO.finish!(p)  # would set false if value were not sticky
    @test GO.predicate_value(p) == true
end

@testset "is_intersection" begin
    @test GO.is_intersection(GO.LOC_INTERIOR, GO.LOC_INTERIOR)
    @test GO.is_intersection(GO.LOC_BOUNDARY, GO.LOC_INTERIOR)
    @test !GO.is_intersection(GO.LOC_EXTERIOR, GO.LOC_INTERIOR)
    @test !GO.is_intersection(GO.LOC_INTERIOR, GO.LOC_EXTERIOR)
    @test !GO.is_intersection(GO.LOC_EXTERIOR, GO.LOC_EXTERIOR)
end

@testset "IMPredicate core state" begin
    # No named IM kinds exist until Task 4; exercise the shared core
    # through a minimal local kind.
    struct _TestIMKind end
    # never determined early; value is whether I/I intersects
    GO.is_determined(p::GO.IMPredicate{_TestIMKind}) = false
    GO.value_im(p::GO.IMPredicate{_TestIMKind}) =
        GO.is_intersects_entry(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR)

    p = GO.IMPredicate(_TestIMKind())
    # JTS IMPredicate constructor presets E/E to dim 2 (Dimension.A)
    @test GO.get_dimension(p, GO.LOC_EXTERIOR, GO.LOC_EXTERIOR) == GO.DIM_A
    # all other entries start FALSE (JTS IntersectionMatrix init)
    @test GO.get_dimension(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR) == GO.DIM_FALSE
    @test !GO.is_known(p)

    GO.init_dims!(p, 2, 1)
    @test p.dimA == GO.DIM_A && p.dimB == GO.DIM_L

    # only an increased dimension is recorded
    @test GO.is_dim_changed(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_L)
    @test GO.get_dimension(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR) == GO.DIM_L
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_P)
    @test GO.get_dimension(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR) == GO.DIM_L
    @test GO.is_dimension_entry(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_L)
    @test GO.is_intersects_entry(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR)
    @test !GO.is_intersects_entry(p, GO.LOC_BOUNDARY, GO.LOC_BOUNDARY)

    # intersects_exterior_of
    @test !GO.intersects_exterior_of(p, true)
    @test !GO.intersects_exterior_of(p, false)
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L)
    @test GO.intersects_exterior_of(p, false)
    @test !GO.intersects_exterior_of(p, true)

    # finish! resolves the value from the matrix state
    @test !GO.is_known(p)
    GO.finish!(p)
    @test GO.is_known(p) && GO.predicate_value(p) == true

    # is_dims_compatible_with_covers (JTS IMPredicate static helper)
    @test GO.is_dims_compatible_with_covers(GO.DIM_P, GO.DIM_L)  # points coveredBy zero-length lines
    @test GO.is_dims_compatible_with_covers(GO.DIM_A, GO.DIM_L)
    @test GO.is_dims_compatible_with_covers(GO.DIM_L, GO.DIM_L)
    @test !GO.is_dims_compatible_with_covers(GO.DIM_L, GO.DIM_A)
end
