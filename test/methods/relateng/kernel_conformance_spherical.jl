# Spherical RelateKernel conformance suite (design layer contract, Task 9).
# Standalone counterpart to kernel_conformance.jl. Inputs are exactly-
# representable integer xyz vectors in *general position*: every predicate is
# a sign of an integer determinant, so the exact answer is unambiguous and the
# suite is deterministic. Points are NOT unit length — sign predicates do not
# require it; where unit length matters (arc membership) we use exact-on-grid
# configurations.

using Test
import GeometryOps as GO
import GeometryOps: Spherical, True, False
import GeometryOps.UnitSpherical: UnitSphericalPoint
import GeoInterface as GI
using Random
using LinearAlgebra: ⋅, cross

const USPt = UnitSphericalPoint{Float64}
_sgn(x) = Int(sign(x))
_usp(x, y, z) = UnitSphericalPoint(Float64(x), Float64(y), Float64(z))

function kernel_conformance_suite_spherical(m; exact)
    rng = Random.MersenneTwister(0x5e1a7e)
    # random integer-component direction (a point on the sphere only up to
    # scaling; sign predicates are scale-invariant in each argument)
    rpt() = _usp(rand(rng, -8:8), rand(rng, -8:8), rand(rng, -8:8))
    function nonzero()
        p = rpt()
        while iszero(GI.x(p)) && iszero(GI.y(p)) && iszero(GI.z(p))
            p = rpt()
        end
        return p
    end

    @testset "rk_orient: antisymmetry / cyclic invariance / degeneracy" begin
        for _ in 1:500
            a, b, c = nonzero(), nonzero(), nonzero()
            o = _sgn(GO.rk_orient(m, a, b, c; exact))
            @test o == -_sgn(GO.rk_orient(m, b, a, c; exact))
            @test o == -_sgn(GO.rk_orient(m, a, c, b; exact))
            @test o == _sgn(GO.rk_orient(m, b, c, a; exact))
            @test o == _sgn(GO.rk_orient(m, c, a, b; exact))
            @test GO.rk_orient(m, a, a, b; exact) == 0
            @test GO.rk_orient(m, a, b, a; exact) == 0
            @test GO.rk_orient(m, a, b, b; exact) == 0
        end
    end
    # testsets added task-by-task below
end

@testset "Kernel conformance: Spherical" begin
    @testset "exact = $E" for E in (True(), False())
        kernel_conformance_suite_spherical(Spherical(); exact = E)
    end
end
