module Predicates

    using ExactPredicates
    using ExactPredicates.Codegen
    using ExactPredicates.Codegen: group!, coord, @genpredicate

    using ExactPredicates: ext, inp

    # This is the implementation of r_cross_s from 2 points in the `_intersection_points` file.
    # 0 == parallel
    function isparallel(a1, a2, b1, b2)
        r = a2 - a1
        s = b2 - b1
        isparallel(r, s)
    end

    @genpredicate function isparallel(r :: 2, s :: 2)
        group!(r...)
        group!(s...)
        ext(r, s)
    end
    # This is the implementation of `iscollinear` from `intersection_points`.
    # 0 == parallel.
    function iscollinear(a1, a2, b1, b2)
        Δqp = b1 - a1
        s = b2 - b1
        iscollinear(Δqp, s)
    end

    @genpredicate function iscollinear(Δqp :: 2, s :: 2)
        group!(Δqp...)
        group!(s...)
        ext(s, Δqp)
    end 

    function istouching(a1, a2, b1, b2)
        pqa = orient(p, q, a)
        pqb = orient(p, q, b)
        touch = if pqa == 0 || pqb == 0
            0
        elseif opposite_signs(pqa, pqb)
            1
        else
            -1
        end
        return touch
    end
end

import .Predicates

# Predicates.r_cross_s(...)

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