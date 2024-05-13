module Predicates
    using ExactPredicates, ExactPredicates.Codegen
    import ExactPredicates: ext
    import ExactPredicates.Codegen: group!, coord, @genpredicate
    import GeometryOps: _False, _True, _booltype
    import GeoInterface as GI

    orient(args...) = ExactPredicates.orient(args...)

    sameside(args...) = ExactPredicates.sameside(args...)

    cross(a, b; exact = _False()) = cross(_booltype(exact), a, b)

    function cross(::_False, a, b)
        c_t1 = GI.x(a) * GI.y(b)
        c_t2 = GI.y(a) * GI.x(b)
        c_val = if c_t1 ≈ c_t2
            0
        else
            c = c_t1 - c_t2
            c > 0 ? 1 : -1
        end
        return c_val
    end

    cross(::_True, a, b) = cross_exact(a, b)

    @genpredicate function cross_exact(a :: 2, b :: 2)
        group!(a...)
        group!(b...)
        ext(a, b)
    end 
end

import .Predicates

#=
# If we want to inject adaptivity, we would do:
using MultiFloats
function isparallel_adaptive(a1, a2, b1, b2)
    r = Float64x2.(a2) - Float64x2.(a1)
    s = Float64x2.(b2) - Float64x2.(b1)
    return sign(Predicates.ext(r, s))
end
function isparallel(a1, a2, b1, b2)
    # try Predicates.isparallel_naive(a1, a2, b1, b2)
    # check the error bound there
    # then try isparallel_adaptive(a1, a2, b1, b2)
    # then try Predicates.isparallel_slow
    # then try isparallel
end
=#