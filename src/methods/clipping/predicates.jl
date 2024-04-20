module Predicates

using ExactPredicates
using ExactPredicates.Codegen
using ExactPredicates.Codegen: group!, coord, @genpredicate

using ExactPredicates: ext, inp

# This is the implementation of r_cross_s from 2 points in the `_intersection_points` file.
@genpredicate function isparallel(a1 :: 2, a2 :: 2, b1 :: 2, b2 :: 2)
    r = a2 - a1
    s = b2 - b1
    group!(r...)
    group!(s...)
    ext(r, s)
end

@genpredicate function iscollinear(a1::2, a2::2, b1::2, b2::2)
    Δqp = b1 - a1
    s = b2 - b1
    group!(Δqp...)
    group!(s...)
    ext(s, Δqp)
end

import .Predicates

# Predicates.r_cross_s(...)

#=
using MultiFloats
function isparallel_adaptive(a1, a2, b1, b2)
    r = Float64x2.(a2) - Float64x2.(a1)
    s = Float64x2.(b2) - Float64x2.(b1)
    return sign(Predicates.ext(r, s))
end
=#