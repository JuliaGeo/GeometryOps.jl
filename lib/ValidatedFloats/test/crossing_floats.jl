module TestCrossingFloats

using Test
using Random
using LinearAlgebra
using StaticArrays
import ValidatedFloats
const CF = ValidatedFloats.CrossingFloats
const R = Rational{BigInt}

encloses(x, exact) = !CF.isbounded(x) ||
    abs(R(x.hi) + R(x.lo) - exact) <= R(CF.radius(x))

struct StaticPoint{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

@testset "restricted scalar and fused vectors" begin
    a = CF.CrossingFloat.((0.8, 0.6, 0.0))
    b = CF.CrossingFloat.((0.8, 0.0, 0.6))
    c = cross(a, b)
    exact = (R(0.6) * R(0.6), -R(0.8) * R(0.6), -R(0.6) * R(0.8))
    @test all(encloses(c[i], exact[i]) for i in 1:3)
    @test encloses(dot(a, b), sum(R(a[i].hi) * R(b[i].hi) for i in 1:3))
end

@testset "power-of-two normalized division" begin
    for (a, b) in ((0x1.4p-200, 0x1.8p-300),
                   (-0x1.cp250, 0x1.2p100),
                   (0x1.1p-250, -0x1.ap-100),
                   (0x1.8p300, 0x1.4p200))
        q = @inferred CF.CrossingFloat(a) / CF.CrossingFloat(b)
        @test CF.isbounded(q)
        @test encloses(q, R(a) / R(b))
    end

    # Cancellation produces a nonzero low limb before division; normalization
    # must preserve its contribution to the enclosure.
    x = CF.CrossingFloat(1.0) + CF.CrossingFloat(0x1p-52) - CF.CrossingFloat(1.0)
    y = CF.CrossingFloat(0x1.8p-200)
    @test encloses(@inferred(x / y), R(0x1p-52) / R(0x1.8p-200))

    # Propagated input radii remain enclosed after a large denominator shift.
    xr = CF.raw(0x1.4p20, 0x1p-34, 0x1p-40)
    yr = CF.raw(0x1.8p-180, 0x1p-234, 0x1p-240)
    q = xr / yr
    @test CF.isbounded(q)
    for xq in (R(xr.hi) + R(xr.lo) - R(xr.rad),
               R(xr.hi) + R(xr.lo) + R(xr.rad)),
        yq in (R(yr.hi) + R(yr.lo) - R(yr.rad),
               R(yr.hi) + R(yr.lo) + R(yr.rad))
        @test encloses(q, xq / yq)
    end

    @test !CF.isbounded(CF.CrossingFloat(1.0) / CF.CrossingFloat(0.0))
    @test !CF.isbounded(CF.CrossingFloat(0x1p400) / CF.CrossingFloat(0x1p-400))
    uncertain = CF.raw(0x1p-200, 0.0, 0.2 * 0x1p-200)
    @test !CF.isbounded(CF.CrossingFloat(1.0) / uncertain)

    rng = MersenneTwister(0xd1a1de)
    for _ in 1:1000
        ey = rand(rng, -300:300)
        eq = rand(rng, -80:80)
        a = ldexp(1.0 + rand(rng), ey + eq)
        b = copysign(ldexp(1.0 + rand(rng), ey), rand(rng, Bool) ? 1.0 : -1.0)
        q = CF.CrossingFloat(a) / CF.CrossingFloat(b)
        @test CF.isbounded(q)
        @test encloses(q, R(a) / R(b))
    end
end

@testset "StaticArrays LinearAlgebra interface" begin
    a = SVector(CF.CrossingFloat.((0.8, 0.6, 0.0)))
    b = SVector(CF.CrossingFloat.((0.8, 0.0, 0.6)))
    c = @inferred cross(a, b)
    d = @inferred dot(a, b)
    s = @inferred a + b
    @test c isa SVector{3,CF.CrossingFloat}
    @test d isa CF.CrossingFloat
    @test s isa SVector{3,CF.CrossingFloat}
    @test c == SVector(cross(Tuple(a), Tuple(b)))
    @test d == dot(Tuple(a), Tuple(b))
    @test which(cross, (typeof(a), typeof(b))).module === CF
    @test which(dot, (typeof(a), typeof(b))).module === CF
    exact_cross = (R(0.6) * R(0.6), -R(0.8) * R(0.6), -R(0.6) * R(0.8))
    @test all(encloses(c[i], exact_cross[i]) for i in 1:3)
    @test encloses(d, sum(R(a[i].hi) * R(b[i].hi) for i in 1:3))
end

@testset "generic FieldVector LinearAlgebra interface" begin
    a = StaticPoint(CF.CrossingFloat(0.8), CF.CrossingFloat(0.6), CF.CrossingFloat(0.0))
    b = StaticPoint(CF.CrossingFloat(0.8), CF.CrossingFloat(0.0), CF.CrossingFloat(0.6))
    c = @inferred cross(a, b)
    d = @inferred dot(a, b)
    @test c isa SVector{3,CF.CrossingFloat}
    @test d isa CF.CrossingFloat
    @test which(cross, (typeof(a), typeof(b))).module === CF
    @test which(dot, (typeof(a), typeof(b))).module === CF
    exact_cross = (R(0.6) * R(0.6), -R(0.8) * R(0.6), -R(0.6) * R(0.8))
    @test all(encloses(c[i], exact_cross[i]) for i in 1:3)
    @test encloses(d, sum(R(a[i].hi) * R(b[i].hi) for i in 1:3))
end

@testset "tuple LinearAlgebra interface" begin
    a = CF.CrossingFloat.((0.8, 0.6, 0.0))
    b = CF.CrossingFloat.((0.8, 0.0, 0.6))
    c = @inferred cross(a, b)
    d = @inferred dot(a, b)
    @test c isa NTuple{3,CF.CrossingFloat}
    @test d isa CF.CrossingFloat
    exact_cross = (R(0.6) * R(0.6), -R(0.8) * R(0.6), -R(0.6) * R(0.8))
    @test all(encloses(c[i], exact_cross[i]) for i in 1:3)
    @test encloses(d, sum(R(a[i].hi) * R(b[i].hi) for i in 1:3))
    @test which(cross, (typeof(a), typeof(b))).module === CF
    @test which(dot, (typeof(a), typeof(b))).module === CF
    for n in (2, 4)
        x = ntuple(i -> CF.CrossingFloat(Float64(i)), n)
        @test encloses(@inferred(dot(x, x)), sum(R(i)^2 for i in 1:n))
    end
end

@testset "Exact3 admission and coordinate planes" begin
    @test all(name -> name ∉ names(CF), (:Exact3, :exact3, :cross3, :dot3, :sum3))
    ex = CF.exact3((1.0, 0.0, 0.0))
    ey = CF.exact3((0.0, 1.0, 0.0))
    @test ex isa CF.Exact3
    @test CF.exact3((Inf, 0.0, 0.0)) === nothing
    @test CF.exact3((0x1p-500, 0.0, 0.0)) === nothing
    @test_throws MethodError CF.Exact3((1.0, 2.0, 3.0))
    c = CF.cross3(ex, ey)
    @test CF.certify.(Ref(Float64), c) == (0.0, 0.0, 1.0)
    @test CF.certify.(Ref(Float64), CF.sum3(ex, ey)) == (1.0, 1.0, 0.0)
end

@testset "fused and Exact3 exact enclosure" begin
    rng = MersenneTwister(0x51a7)
    for _ in 1:2000
        a = ntuple(_ -> randn(rng), 3)
        b = ntuple(_ -> randn(rng), 3)
        ea, eb = CF.exact3(a), CF.exact3(b)
        @test ea isa CF.Exact3
        @test eb isa CF.Exact3
        ce = CF.cross3(ea, eb)
        cs = CF.cross3(CF.CrossingFloat.(a), CF.CrossingFloat.(b))
        exact = (R(a[2]) * R(b[3]) - R(a[3]) * R(b[2]),
                 R(a[3]) * R(b[1]) - R(a[1]) * R(b[3]),
                 R(a[1]) * R(b[2]) - R(a[2]) * R(b[1]))
        for i in 1:3
            @test encloses(ce[i], exact[i])
            @test encloses(cs[i], exact[i])
        end
    end
end

end
