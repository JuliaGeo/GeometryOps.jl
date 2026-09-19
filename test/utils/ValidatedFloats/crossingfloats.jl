using Test, Random
using .CrossingFloats
const CF = CrossingFloats
const CQ = Rational{BigInt}

function cfcontains(x, q::CQ)
    !CF.isbounded(x) || abs(CQ(x.hi)+CQ(x.lo)-q) <= CQ(x.rad)
end
function cfinterval(x)
    c = CQ(x.hi)+CQ(x.lo)
    return c-CQ(x.rad),c+CQ(x.rad)
end

@testset "restricted crossing arithmetic" begin
    rng = MersenneTwister(511507)
    @test CF.certify(Float64,CF.CrossingFloat(0.0)) == 0.0
    @test !CF.isbounded(CF.CrossingFloat(Inf))
    @test !CF.isbounded(CF.CrossingFloat(0x1p-500))
    @test !CF.isbounded(CF.CrossingFloat(0x1p300)*CF.CrossingFloat(0x1p300))
    for _ in 1:2000
        a,b,c = ntuple(_ -> ldexp(randn(rng),rand(rng,-180:180)),3)
        x,y,z = CF.CrossingFloat.((a,b,c))
        ra,rb,rc = CQ.((a,b,c))
        @test cfcontains(x+y,ra+rb)
        @test cfcontains(x*y,ra*rb)
        @test cfcontains((x+y)*z,(ra+rb)*rc)
        @test cfcontains(x*y-x*z,ra*rb-ra*rc)
        d = 1.0+3rand(rng)
        @test cfcontains((x+y)/CF.CrossingFloat(d),(ra+rb)/CQ(d))
    end
    for _ in 1:1000
        # Exercise propagation of nonzero input radii, including their product.
        ah,bh = 1+2rand(rng),1+2rand(rng)
        al,bl = (rand(rng)-0.5)*eps(ah),(rand(rng)-0.5)*eps(bh)
        ar,br = ldexp(ah,rand(rng,-100:-20)),ldexp(bh,rand(rng,-100:-20))
        a,b = CF.raw(ah,al,ar),CF.raw(bh,bl,br)
        for qa in cfinterval(a), qb in cfinterval(b)
            @test cfcontains(a+b,qa+qb)
            @test cfcontains(a*b,qa*qb)
            @test cfcontains(a/b,qa/qb)
        end
        s = sqrt(a)
        @test CF.isbounded(s)
        slo,shi = cfinterval(s)
        alo,ahi = cfinterval(a)
        @test slo >= 0 && slo*slo <= alo
        @test shi*shi >= ahi
    end
    for e in (-399,-100,0,100,399), sg in (-1.0,1.0)
        f = sg*ldexp(1.0,e)
        x = CF.CrossingFloat(f)
        @test CF.certify(Float64,x) == f
        # A radius spanning an adjacent rounding bin must be rejected.
        @test CF.certify(Float64,CF.raw(f,0.0,eps(f))) === nothing
    end
end
