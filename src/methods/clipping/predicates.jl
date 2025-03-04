module Predicates
    using ExactPredicates, ExactPredicates.Codegen
    import ExactPredicates: ext
    import ExactPredicates.Codegen: group!, @genpredicate
    import GeometryOps: False, True, booltype, _tuple_point
    import GeoInterface as GI
    import AdaptivePredicates

    #= Determine the orientation of c with regards to the oriented segment (a, b).
    Return 1 if c is to the left of (a, b).
    Return -1 if c is to the right of (a, b).
    Return 0 if c is on (a, b) or if a == b. =#
    orient(a, b, c; exact) = _orient(booltype(exact), _tuple_point(a, Float64), _tuple_point(b, Float64), _tuple_point(c, Float64))
    
    # If `exact` is `true`, use `ExactPredicates` to calculate the orientation.
    _orient(::True, a, b, c) = AdaptivePredicates.orient2p(_tuple_point(a, Float64), _tuple_point(b, Float64), _tuple_point(c, Float64))
    # _orient(::True, a, b, c) = ExactPredicates.orient(_tuple_point(a, Float64), _tuple_point(b, Float64), _tuple_point(c, Float64))
    # If `exact` is `false`, calculate the orientation without using `ExactPredicates`.
    function _orient(exact::False, a, b, c)
        a = a .- c
        b = b .- c
        return _cross(exact, a, b)
    end

    #= Determine the sign of the cross product of a and b.
    Return 1 if the cross product is positive.
    Return -1 if the cross product is negative.
    Return 0 if the cross product is 0. =#
    cross(a, b; exact) = _cross(booltype(exact), a, b)

    #= If `exact` is `true`, use exact cross product calculation created using
    `ExactPredicates`generated predicate. Note that as of now `ExactPredicates` requires
    Float64 so we must convert points a and b. =#
    _cross(::True, a, b) = _cross_exact(_tuple_point(a, Float64), _tuple_point(b, Float64))

    # Exact cross product calculation using `ExactPredicates`.
    @genpredicate function _cross_exact(a :: 2, b :: 2)
        group!(a...)
        group!(b...)
        ext(a, b)
    end 

    # If `exact` is `false`, calculate the cross product without using `ExactPredicates`.
    function _cross(::False, a, b)
        c_t1 = GI.x(a) * GI.y(b)
        c_t2 = GI.y(a) * GI.x(b)
        c_val = if isapprox(c_t1, c_t2)
            0
        else
            sign(c_t1 - c_t2)
        end
        return c_val
    end
    
end

import .Predicates

#=
# If we want to inject adaptivity, we would do something like:
function cross(a, b, c)
    # try Predicates._cross_naive(a, b, c)
    # check the error bound there
    # then try Predicates._cross_adaptive(a, b, c)
    # then try Predicates._cross_exact
end
=#