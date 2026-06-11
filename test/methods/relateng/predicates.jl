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

# --- Named IM predicates, pattern matcher, matrix predicate (Task 4) ---
# Port of JTS RelatePredicateTest.java, plus requirement-flag-table
# assertions read directly from each RelatePredicate.java inner class,
# and IMPatternMatcher / RelateMatrixPredicate tests.

# JTS RelatePredicateTest fixture matrices ('.' separators are cosmetic)
const A_EXT_B_INT = "***.***.1**"
const A_INT_B_INT = "1**.***.***"

# JTS RelatePredicateTest.applyIM
function apply_im!(pred, im_in::String)
    im = replace(im_in, "." => "")
    for i in 0:8
        locA = i ÷ 3
        locB = i - 3 * locA
        entry = im[i + 1]
        if entry in ('0', '1', '2')
            GO.update_dim!(pred, locA, locB, GO.dim_code(entry))
        end
    end
    return pred
end

# JTS RelatePredicateTest.checkPredicate / checkPred
function check_predicate(pred, im::String, expected::Bool)
    apply_im!(pred, im)
    GO.finish!(pred)
    @test GO.is_known(pred)
    @test GO.predicate_value(pred) == expected
end

# JTS RelatePredicateTest.checkPredicatePartial (value known before finish!)
function check_predicate_partial(pred, im::String, expected::Bool)
    apply_im!(pred, im)
    @test GO.is_known(pred)  # "predicate value is not known"
    GO.finish!(pred)
    @test GO.predicate_value(pred) == expected
end

@testset "RelatePredicateTest (JTS port)" begin
    @testset "testIntersects" begin
        check_predicate(GO.pred_intersects(), A_INT_B_INT, true)
    end
    @testset "testDisjoint" begin
        check_predicate(GO.pred_intersects(), A_EXT_B_INT, false)
        check_predicate(GO.pred_disjoint(), A_EXT_B_INT, true)
    end
    @testset "testCovers" begin
        check_predicate(GO.pred_covers(), A_INT_B_INT, true)
        check_predicate(GO.pred_covers(), A_EXT_B_INT, false)
    end
    @testset "testCoversFast" begin
        check_predicate_partial(GO.pred_covers(), A_EXT_B_INT, false)
    end
    @testset "testMatch" begin
        check_predicate(GO.pred_matches("1***T*0**"), "1**.*2*.0**", true)
    end
end

@testset "predicate names" begin
    @test GO.predicate_name(GO.pred_contains()) == "contains"
    @test GO.predicate_name(GO.pred_within()) == "within"
    @test GO.predicate_name(GO.pred_covers()) == "covers"
    @test GO.predicate_name(GO.pred_coveredby()) == "coveredBy"
    @test GO.predicate_name(GO.pred_crosses()) == "crosses"
    @test GO.predicate_name(GO.pred_equalstopo()) == "equals"
    @test GO.predicate_name(GO.pred_overlaps()) == "overlaps"
    @test GO.predicate_name(GO.pred_touches()) == "touches"
end

@testset "requirement flag table (per RelatePredicate.java)" begin
    # columns: self_noding, interaction, covers(A), covers(B), ext_check(A), ext_check(B)
    flag_table = [
        (GO.pred_contains(),   true, true,  true,  false, false, true),
        (GO.pred_within(),     true, true,  false, true,  true,  false),
        (GO.pred_covers(),     true, true,  true,  false, false, true),
        (GO.pred_coveredby(),  true, true,  false, true,  true,  false),
        (GO.pred_crosses(),    true, true,  false, false, true,  true),
        (GO.pred_equalstopo(), true, false, false, false, true,  true),
        (GO.pred_overlaps(),   true, true,  false, false, true,  true),
        (GO.pred_touches(),    true, true,  false, false, true,  true),
    ]
    for (p, sn, ia, cov_a, cov_b, ext_a, ext_b) in flag_table
        T = typeof(p)
        @testset "$(GO.predicate_name(p))" begin
            @test GO.require_self_noding(T) == sn
            @test GO.require_interaction(T) == ia
            @test GO.require_covers(T, true) == cov_a
            @test GO.require_covers(T, false) == cov_b
            @test GO.require_exterior_check(T, true) == ext_a
            @test GO.require_exterior_check(T, false) == ext_b
            # instance-level forwarding agrees with the type-level flags
            @test GO.require_self_noding(p) == sn
            @test GO.require_interaction(p) == ia
            @test GO.require_covers(p, true) == cov_a && GO.require_covers(p, false) == cov_b
            @test GO.require_exterior_check(p, true) == ext_a && GO.require_exterior_check(p, false) == ext_b
        end
    end
end

@testset "contains" begin
    # dims incompatible with covers determine false at init_dims!
    p = GO.pred_contains()
    GO.init_dims!(p, GO.DIM_L, GO.DIM_A)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # B in interior of A
    check_predicate(GO.pred_contains(), A_INT_B_INT, true)
    # B intersecting exterior of A determines false early
    check_predicate_partial(GO.pred_contains(), A_EXT_B_INT, false)
    # envelope of A must cover envelope of B
    q = GO.pred_contains()
    GO.init_bounds!(q, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(0.0, 2.0), Y=(0.0, 1.0)))
    @test GO.is_known(q) && GO.predicate_value(q) == false
    r = GO.pred_contains()
    GO.init_bounds!(r, Extents.Extent(X=(0.0, 3.0), Y=(0.0, 3.0)),
                       Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0)))
    @test !GO.is_known(r)
end

@testset "within" begin
    p = GO.pred_within()
    GO.init_dims!(p, GO.DIM_A, GO.DIM_L)   # dimB must be able to cover dimA
    @test GO.is_known(p) && GO.predicate_value(p) == false
    check_predicate(GO.pred_within(), A_INT_B_INT, true)
    # A intersecting exterior of B determines false early
    check_predicate_partial(GO.pred_within(), "1*1.***.***", false)
    # envelope of B must cover envelope of A
    q = GO.pred_within()
    GO.init_bounds!(q, Extents.Extent(X=(0.0, 2.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)))
    @test GO.is_known(q) && GO.predicate_value(q) == false
end

@testset "covers and coveredBy boundary contact" begin
    # covers is true for boundary-only contact, where contains is false
    check_predicate(GO.pred_covers(), "***.1**.***", true)
    check_predicate(GO.pred_contains(), "***.1**.***", false)
    check_predicate(GO.pred_coveredby(), "**F.1*F.***", true)
    check_predicate_partial(GO.pred_coveredby(), "1*1.***.***", false)
    check_predicate(GO.pred_coveredby(), A_INT_B_INT, true)
end

@testset "crosses" begin
    # P/P and A/A are determined false at init_dims!
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_P, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_A, GO.DIM_A)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # dimA < dimB: determined once I/I and I/E intersect
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_L, GO.DIM_A)
    apply_im!(p, "1*1.***.***")
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # dimA > dimB: determined once I/I and E/I intersect
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_A, GO.DIM_L)
    apply_im!(p, "1**.***.1**")
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # L/L: dim-0 interior intersection is a crossing, found at finish!
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_L, GO.DIM_L)
    apply_im!(p, "0**.***.***")
    @test !GO.is_known(p)
    GO.finish!(p)
    @test GO.predicate_value(p) == true
    # L/L: dim-1 (collinear) interior intersection is not a crossing
    p = GO.pred_crosses(); GO.init_dims!(p, GO.DIM_L, GO.DIM_L)
    apply_im!(p, "1**.***.***")
    @test GO.is_known(p) && GO.predicate_value(p) == false
end

@testset "equalsTopo" begin
    p = GO.pred_equalstopo(); GO.init_dims!(p, GO.DIM_L, GO.DIM_L)
    apply_im!(p, A_INT_B_INT)
    GO.finish!(p)
    @test GO.predicate_value(p) == true
    # different dims are never topo-equal
    p = GO.pred_equalstopo(); GO.init_dims!(p, GO.DIM_L, GO.DIM_A)
    apply_im!(p, A_INT_B_INT)
    GO.finish!(p)
    @test GO.predicate_value(p) == false
    # any exterior intersection determines false early
    p = GO.pred_equalstopo(); GO.init_dims!(p, GO.DIM_A, GO.DIM_A)
    apply_im!(p, "2*1.***.***")
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # EMPTY = EMPTY (null bounds)
    p = GO.pred_equalstopo()
    GO.init_bounds!(p, nothing, nothing)
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # unequal bounds determine false
    p = GO.pred_equalstopo()
    GO.init_bounds!(p, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(0.0, 2.0), Y=(0.0, 1.0)))
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # equal bounds leave the value unknown
    p = GO.pred_equalstopo()
    GO.init_bounds!(p, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)))
    @test !GO.is_known(p)
end

@testset "overlaps" begin
    # different dims are determined false at init_dims!
    p = GO.pred_overlaps(); GO.init_dims!(p, GO.DIM_L, GO.DIM_A)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # A/A: determined once I/I, I/E and E/I all intersect
    p = GO.pred_overlaps(); GO.init_dims!(p, GO.DIM_A, GO.DIM_A)
    apply_im!(p, "2*2.***.2**")
    @test GO.is_known(p) && GO.predicate_value(p) == true
    # L/L: requires a dim-1 interior/interior intersection
    p = GO.pred_overlaps(); GO.init_dims!(p, GO.DIM_L, GO.DIM_L)
    apply_im!(p, "0*1.***.1**")
    GO.finish!(p)
    @test GO.predicate_value(p) == false
    p = GO.pred_overlaps(); GO.init_dims!(p, GO.DIM_L, GO.DIM_L)
    apply_im!(p, "1*1.***.1**")
    @test GO.is_known(p) && GO.predicate_value(p) == true
end

@testset "touches" begin
    # Points have only interiors, so cannot touch
    p = GO.pred_touches(); GO.init_dims!(p, GO.DIM_P, GO.DIM_P)
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # boundary-only contact touches
    p = GO.pred_touches(); GO.init_dims!(p, GO.DIM_A, GO.DIM_A)
    apply_im!(p, "***.*1*.***")
    GO.finish!(p)
    @test GO.predicate_value(p) == true
    # interior intersection determines false early
    p = GO.pred_touches(); GO.init_dims!(p, GO.DIM_A, GO.DIM_A)
    apply_im!(p, "2**.***.***")
    @test GO.is_known(p) && GO.predicate_value(p) == false
    # dim order is symmetric (JTS isTouches transposes dims, not the matrix)
    p = GO.pred_touches(); GO.init_dims!(p, GO.DIM_A, GO.DIM_L)
    apply_im!(p, "***.*1*.***")
    GO.finish!(p)
    @test GO.predicate_value(p) == true
end

@testset "DE9IM relate queries (IntersectionMatrix port)" begin
    contains_im = GO.DE9IM("212FF1FF2")
    within_im   = GO.DE9IM("2FF1FF212")
    equals_im   = GO.DE9IM("2FF1FFFF2")
    @test GO.is_contains(contains_im) && !GO.is_contains(within_im)
    @test GO.is_within(within_im) && !GO.is_within(contains_im)
    @test GO.is_covers(contains_im) && !GO.is_covers(within_im)
    @test GO.is_coveredby(within_im) && !GO.is_coveredby(contains_im)
    # covers via boundary contact only, where contains fails
    boundary_im = GO.DE9IM("FF2F01FF2")
    @test GO.is_covers(boundary_im) && !GO.is_contains(boundary_im)
    # equals requires equal dims and no exterior intersections
    @test GO.is_equals(equals_im, GO.DIM_A, GO.DIM_A)
    @test !GO.is_equals(equals_im, GO.DIM_A, GO.DIM_L)
    @test !GO.is_equals(contains_im, GO.DIM_A, GO.DIM_A)
    # crosses
    la_cross = GO.DE9IM("101FF0212")
    @test GO.is_crosses(la_cross, GO.DIM_L, GO.DIM_A)
    @test GO.is_crosses(la_cross, GO.DIM_A, GO.DIM_L)
    @test !GO.is_crosses(la_cross, GO.DIM_L, GO.DIM_L)   # L/L needs I/I == 0
    @test !GO.is_crosses(la_cross, GO.DIM_A, GO.DIM_A)   # A/A never crosses
    @test GO.is_crosses(GO.DE9IM("0F1FF0102"), GO.DIM_L, GO.DIM_L)
    # overlaps
    @test GO.is_overlaps(GO.DE9IM("212101212"), GO.DIM_A, GO.DIM_A)
    @test !GO.is_overlaps(GO.DE9IM("212101212"), GO.DIM_L, GO.DIM_L)
    @test GO.is_overlaps(GO.DE9IM("1F1FF0102"), GO.DIM_L, GO.DIM_L)
    @test !GO.is_overlaps(GO.DE9IM("212101212"), GO.DIM_L, GO.DIM_A)
    # touches
    touch_im = GO.DE9IM("FF2F11212")
    @test GO.is_touches(touch_im, GO.DIM_A, GO.DIM_A)
    @test GO.is_touches(touch_im, GO.DIM_A, GO.DIM_L)  # transposed dims
    @test !GO.is_touches(touch_im, GO.DIM_P, GO.DIM_P) # points cannot touch
    @test !GO.is_touches(GO.DE9IM("212101212"), GO.DIM_A, GO.DIM_A)
end

@testset "IMPatternMatcher" begin
    p = GO.IMPatternMatcher("T*F**FFF*")
    @test GO.predicate_name(p) == "IMPattern"
    @test !GO.is_known(p)
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    @test !GO.is_known(p)
    # an entry violating the pattern mask determines false immediately
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_L)  # pattern 'F' at I/E
    @test GO.is_known(p) && GO.predicate_value(p) == false

    # a fully matching evaluation
    q = GO.IMPatternMatcher("T*F**FFF*")
    GO.update_dim!(q, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    GO.update_dim!(q, GO.LOC_BOUNDARY, GO.LOC_INTERIOR, GO.DIM_L)
    GO.finish!(q)
    @test GO.is_known(q) && GO.predicate_value(q) == true

    # pred_matches is the JTS RelatePredicate.matches factory
    @test GO.pred_matches("T*F**FFF*") isa GO.IMPatternMatcher

    # require_interaction is an instance method computed from the pattern
    @test GO.require_interaction(GO.IMPatternMatcher("T*F**FFF*")) == true
    @test GO.require_interaction(GO.IMPatternMatcher("FF*FF****")) == false  # disjoint pattern
    @test GO.require_interaction(GO.IMPatternMatcher("****1****")) == true
    @test GO.require_interaction(GO.IMPatternMatcher("FF*FF*1**")) == false  # only E-row entries

    # init_bounds!: interaction-requiring pattern + disjoint envelopes => false
    r = GO.IMPatternMatcher("T*F**FFF*")
    GO.init_bounds!(r, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(5.0, 6.0), Y=(5.0, 6.0)))
    @test GO.is_known(r) && GO.predicate_value(r) == false
    # non-interaction pattern is not determined by disjoint envelopes
    s = GO.IMPatternMatcher("FF*FF****")
    GO.init_bounds!(s, Extents.Extent(X=(0.0, 1.0), Y=(0.0, 1.0)),
                       Extents.Extent(X=(5.0, 6.0), Y=(5.0, 6.0)))
    @test !GO.is_known(s)
end

@testset "IntersectionMatrixPattern constants" begin
    @test GO.IM_PATTERN_ADJACENT == "F***1****"
    @test GO.IM_PATTERN_CONTAINS_PROPERLY == "T**FF*FF*"
    @test GO.IM_PATTERN_INTERIOR_INTERSECTS == "T********"
end

@testset "RelateMatrixPredicate" begin
    p = GO.RelateMatrixPredicate()
    @test GO.predicate_name(p) == "relateMatrix"
    @test GO.require_interaction(typeof(p)) == false
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_INTERIOR, GO.DIM_A)
    @test !GO.is_known(p)   # never determined early
    GO.update_dim!(p, GO.LOC_INTERIOR, GO.LOC_EXTERIOR, GO.DIM_A)
    @test !GO.is_known(p)
    GO.finish!(p)
    @test !isnothing(GO.result_im(p))
    @test GO.result_im(p)[GO.LOC_INTERIOR, GO.LOC_INTERIOR] == GO.DIM_A
    @test GO.result_im(p)[GO.LOC_INTERIOR, GO.LOC_EXTERIOR] == GO.DIM_A
    # E/E is preset to dim 2 by the IMPredicate constructor
    @test GO.result_im(p)[GO.LOC_EXTERIOR, GO.LOC_EXTERIOR] == GO.DIM_A
end
