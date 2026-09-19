using Test
using GeometryOps.ValidatedFloats: ValidatedFloat, center, radius, isbounded, certify, certified_sign

exactcenter(x) = BigFloat(x.hi) + BigFloat(x.lo)
contains(x, exact) = !isbounded(x) || abs(exactcenter(x) - exact) <= BigFloat(radius(x))

@testset "construction and certification" begin
    for a in (0.0, -0.0, 1.0, -2.0, floatmin(Float64), floatmax(Float64))
        x = ValidatedFloat(a)
        @test contains(x, BigFloat(a))
        @test isequal(certify(Float64, x), a)
    end
    @test !isbounded(ValidatedFloat(Inf))
    @test certify(Float64, ValidatedFloat(Inf)) === nothing
    for hi in (-floatmax(Float64), floatmax(Float64))
        gap = floatmax(Float64) - prevfloat(floatmax(Float64))
        for lo in (-gap, 0.0, gap)
            @test isbounded(ValidatedFloat(hi, lo, 0.0))
        end
        for lo in (-floatmax(Float64), floatmax(Float64), -nextfloat(gap), nextfloat(gap))
            @test_throws ArgumentError ValidatedFloat(hi, lo, 0.0)
        end
    end
end

@testset "rounding boundaries do not falsely certify" begin
    for p in (-100, -1, 0, 1, 100), sign in (-1.0, 1.0)
        f = sign * ldexp(1.0, p)
        toward = sign > 0 ? nextfloat(f) : prevfloat(f)
        halfgap = abs(toward - f) / 2
        midpoint = ValidatedFloat(f, sign * halfgap, 0.0)
        @test certify(Float64, midpoint) === nothing

        # An interval centered on a representable value but reaching a
        # rounding midpoint must also fail, on both sides of zero.
        spanning = ValidatedFloat(f, 0.0, halfgap)
        @test certify(Float64, spanning) === nothing
    end
end

@testset "certified sign" begin
    @test certified_sign(ValidatedFloat(0.0)) == 0
    @test certified_sign(ValidatedFloat(2.0)) == 1
    @test certified_sign(ValidatedFloat(-2.0)) == -1
    @test certified_sign(ValidatedFloat(1.0, 0.0, 1.0)) === nothing
    @test certified_sign(ValidatedFloat(-1.0, 0.0, 1.0)) === nothing
    @test certified_sign(ValidatedFloat(Inf)) === nothing
end

@testset "arithmetic encloses high precision result" begin
    vals = (0.0, 1.0, -1.0, 0.1, -3.25, 0x1p-40, 0x1p40,
            prevfloat(1.0), nextfloat(1.0))
    for a in vals, b in vals
        x, y = ValidatedFloat(a), ValidatedFloat(b)
        @test contains(x + y, BigFloat(a) + BigFloat(b))
        @test contains(x - y, BigFloat(a) - BigFloat(b))
        @test contains(x * y, BigFloat(a) * BigFloat(b))
        b != 0 && @test contains(x / y, BigFloat(a) / BigFloat(b))
    end
    for a in (0.0, 0.1, 1.0, 2.0, 0x1p-40, 0x1p40)
        @test contains(sqrt(ValidatedFloat(a)), sqrt(BigFloat(a)))
    end
end

@testset "cancellation, scaling, and conservative failures" begin
    x = (ValidatedFloat(1.0) + ValidatedFloat(0x1p-53)) - ValidatedFloat(1.0)
    @test contains(x, BigFloat(0x1p-53))
    @test certify(Float64, ValidatedFloat(1.0) - ValidatedFloat(1.0)) == 0.0
    for n in (-900, -100, 0, 100, 900)
        y = ldexp(ValidatedFloat(1.25), n)
        @test contains(y, BigFloat(1.25) * big(2.0)^n)
    end
    @test !isbounded(ValidatedFloat(floatmax(Float64)) * ValidatedFloat(2.0))
    @test !isbounded(ValidatedFloat(floatmin(Float64)) * ValidatedFloat(0.5))
    a = ldexp(nextfloat(1.0), -500)
    z = ValidatedFloat(a) * ValidatedFloat(a)
    @test !isbounded(z)
    @test contains(z, BigFloat(a) * BigFloat(a))
    @test !isbounded(ldexp(ValidatedFloat(floatmin(Float64)), -1))
    @test !isbounded(ldexp(ValidatedFloat(nextfloat(0.0)), -1))
    @test_throws ArgumentError ValidatedFloat(1.0, 1.0, 0.0)
    @test_throws DomainError sqrt(ValidatedFloat(-1.0))
end

@testset "Rational oracle" begin
    for a in (-17:17), b in (1:17)
        x = ValidatedFloat(Float64(a)) / ValidatedFloat(Float64(b))
        @test contains(x, BigFloat(a // b))
    end
end

include("crossing_floats.jl")
include("original_crossingfloats.jl")
