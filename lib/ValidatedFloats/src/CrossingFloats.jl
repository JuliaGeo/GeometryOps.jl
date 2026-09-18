module CrossingFloats

import LinearAlgebra: cross, dot
import StaticArrays: StaticVector, SVector

export CrossingFloat, center, radius, isbounded, certify, certified_sign,
       scale_pow2, PowerOfTwo

# A restricted arithmetic domain keeps the leading-product residual representable.
const LOW = 0x1p-400
const HIGH = 0x1p400
const E = 0x1p-53
const PAD = 1.0 + 0x1p-44
const SHRINK = 1.0 - 0x1p-44
const FLOOR = 0x1p-1000

struct CrossingFloat <: Real
    hi::Float64
    lo::Float64
    rad::Float64
    CrossingFloat(h::Float64, l::Float64, r::Float64, ::Val{:raw}) = new(h,l,r)
end

@inline raw(h,l,r) = CrossingFloat(h,l,r,Val(:raw))
@inline bad() = raw(0.0,0.0,Inf)
@inline finish(::Val{true},h,l,r) = make(h,l,r)
@inline finish(::Val{false},h,l,r) = raw(h,l,r)
@inline admit(x::CrossingFloat) = make(x.hi,x.lo,x.rad)
@inline function admissible(h)
    bits = reinterpret(UInt64,h) & 0x7fffffffffffffff
    return (bits == 0) | ((bits - reinterpret(UInt64,LOW)) <= (reinterpret(UInt64,HIGH)-reinterpret(UInt64,LOW)))
end
@inline make(h,l,r) = admissible(h) & (r <= HIGH) ? raw(h,l,r) : bad()
@inline CrossingFloat(x::Float64) = make(x,0.0,0.0)
@inline center(x::CrossingFloat) = x.hi + x.lo
@inline radius(x::CrossingFloat) = x.rad
@inline isbounded(x::CrossingFloat) = isfinite(x.rad)
@inline pad(r) = (r + FLOOR) * PAD
@inline function twosum(a,b)
    s = a+b; bb = s-a
    return s, (a-(s-bb))+(b-bb)
end
@inline function twoprod(a,b)
    p = a*b
    return p, fma(a,b,-p)
end
@inline Base.:-(x::CrossingFloat) = raw(-x.hi,-x.lo,x.rad)

@inline function add(x::CrossingFloat,y::CrossingFloat,g::Val)
    sh,se = twosum(x.hi,y.hi)
    if (x.lo == 0.0) & (y.lo == 0.0) & (x.rad == 0.0) & (y.rad == 0.0)
        return finish(g,sh,se,0.0)
    end
    t = se+x.lo; t2 = t+y.lo
    rh,rl = twosum(sh,t2)
    b = x.rad+y.rad + E*((abs(se)+abs(x.lo))+(abs(t)+abs(y.lo)))
    return finish(g,rh,rl,pad(b))
end
@inline Base.:+(x::CrossingFloat,y::CrossingFloat) = add(x,y,Val(true))
@inline Base.:-(x::CrossingFloat,y::CrossingFloat) = x+(-y)

@inline function mul(x::CrossingFloat,y::CrossingFloat,g::Val)
    if x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0
        return isbounded(y) ? raw(0.0,0.0,0.0) : bad()
    elseif y.hi == 0.0 && y.lo == 0.0 && y.rad == 0.0
        return isbounded(x) ? raw(0.0,0.0,0.0) : bad()
    end
    ph,pl = twoprod(x.hi,y.hi)
    if (x.lo == 0.0) & (y.lo == 0.0) & (x.rad == 0.0) & (y.rad == 0.0)
        return finish(g,ph,pl,0.0)
    end
    c1 = x.hi*y.lo; c2 = x.lo*y.hi
    c = c1+c2; t = pl+c
    rh,rl = twosum(ph,t)
    er = E*(2*(abs(c1)+abs(c2))+abs(pl)+abs(c)) + abs(x.lo*y.lo)
    b = er + x.rad*(abs(y.hi)+abs(y.lo)+y.rad) + (abs(x.hi)+abs(x.lo))*y.rad
    return finish(g,rh,rl,pad(b))
end

@inline Base.:*(x::CrossingFloat,y::CrossingFloat) = mul(x,y,Val(true))

# The division kernel assumes a normalized denominator. The public method
# reaches that domain with a shared power-of-two scale, which leaves x/y
# unchanged and is exact when `scale_pow2` accepts both operands.
@inline function _div_normalized(x::CrossingFloat,y::CrossingFloat)
    yh = y.hi; yl = y.lo; ey = y.rad
    ((1.0 <= abs(yh) <= 4.0) && abs(yl)+ey <= 0.19abs(yh)) || return bad()
    x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0 && return raw(0.0,0.0,0.0)
    q1 = x.hi/yh
    admissible(q1) || return bad()
    ph,pl = twoprod(q1,yh)
    sh,sl = twosum(x.hi,-ph)
    t1 = sh-pl; t2 = t1+sl; t3 = t2+x.lo
    m = q1*yl; r = t3-m
    q2 = r/yh
    rh,rl = twosum(q1,q2)
    er = E*((abs(sh)+abs(pl))+(abs(t1)+abs(sl))+(abs(t2)+abs(x.lo))+
            abs(m)+(abs(t3)+abs(m)))
    ymin = (abs(yh)-abs(yl))*SHRINK
    ymin2 = (ymin-ey)*SHRINK
    b = E*abs(r)/abs(yh) + abs(r)*abs(yl)/(abs(yh)*ymin) + er/ymin +
        x.rad/ymin2 + (abs(x.hi)+abs(x.lo)+x.rad)*ey/(ymin*ymin2)
    exact = (sh == 0.0) & (sl == 0.0) & (pl == 0.0) & (x.lo == 0.0) &
            (yl == 0.0) & (x.rad == 0.0) & (ey == 0.0)
    return make(rh,rl,exact ? 0.0 : pad(b))
end

@inline function Base.:/(x::CrossingFloat,y::CrossingFloat)
    yh = y.hi
    1.0 <= abs(yh) <= 4.0 && return _div_normalized(x, y)
    isbounded(x) & isbounded(y) || return bad()
    yh == 0.0 && return bad()
    abs(y.lo) + y.rad <= 0.19abs(yh) || return bad()
    p = PowerOfTwo(-exponent(yh))
    xs = scale_pow2(x, p)
    ys = scale_pow2(y, p)
    isbounded(xs) & isbounded(ys) || return bad()
    return _div_normalized(xs, ys)
end

@inline function Base.sqrt(x::CrossingFloat)
    sh,sl,es = x.hi,x.lo,x.rad
    sh == 0.0 && sl == 0.0 && es == 0.0 && return raw(0.0,0.0,0.0)
    ((1.0 <= sh <= 16.0) && abs(sl)+es <= 0.19sh) || return bad()
    s0 = sqrt(sh)
    ph,pl = twoprod(s0,s0)
    r1 = sh-ph; r2 = r1-pl; r = r2+sl
    q = r/(2s0)
    rh,rl = twosum(s0,q)
    g = abs(sl)/(1.8s0) + 2E*s0
    en = g*g/(2s0)
    er = (E*(abs(r1)+abs(r2)+abs(sl)) + E*(abs(r1)+abs(pl)))/(2s0) + E*abs(r)/(2s0)
    b = en+er+es/(0.9s0)
    exact = (r1 == 0.0) & (pl == 0.0) & (sl == 0.0) & (es == 0.0)
    return make(rh,rl,exact ? 0.0 : pad(b))
end

struct PowerOfTwo
    value::Float64
    function PowerOfTwo(n::Integer)
        s = ldexp(1.0,n)
        isfinite(s) && s > 0.0 || throw(ArgumentError("unrepresentable power of two"))
        new(s)
    end
end

@inline function scale_pow2(x::CrossingFloat,p::PowerOfTwo)
    s = p.value
    h,l,r = x.hi*s,x.lo*s,x.rad*s
    ((x.lo == 0.0 || abs(l) >= floatmin(Float64)) &&
     (x.rad == 0.0 || r >= floatmin(Float64)) && (x.hi == 0.0 || h != 0.0)) || return bad()
    return make(h,l,r)
end
@inline Base.ldexp(x::CrossingFloat,n::Integer) = scale_pow2(x,PowerOfTwo(n))

@inline function certified_sign(x::CrossingFloat)
    isbounded(x) || return nothing
    x.hi == 0.0 && x.lo == 0.0 && x.rad == 0.0 && return 0
    abs(x.hi) > (abs(x.lo)+x.rad)*PAD || return nothing
    return signbit(x.hi) ? -1 : 1
end

@inline function certify(::Type{Float64},x::CrossingFloat)
    isbounded(x) || return nothing
    f,rem = twosum(x.hi,x.lo)
    rem == 0.0 && x.rad == 0.0 && return f
    f == 0.0 && return nothing
    gap = abs(f)-prevfloat(abs(f))
    (abs(rem)+x.rad)*PAD < 0.5gap || return nothing
    return f
end

# Products can temporarily exceed the scalar domain; only sums consume them.
# Admitted inputs bound these intermediates below 2^805 and above the EFT grid.
@inline function cross3(a::NTuple{3,CrossingFloat},b::NTuple{3,CrossingFloat})
    g = Val(false)
    return (admit(add(mul(a[2],b[3],g),-mul(a[3],b[2],g),g)),
            admit(add(mul(a[3],b[1],g),-mul(a[1],b[3],g),g)),
            admit(add(mul(a[1],b[2],g),-mul(a[2],b[1],g),g)))
end
@inline function dot3(a::NTuple{3,CrossingFloat},b::NTuple{3,CrossingFloat})
    g = Val(false)
    return admit(add(add(mul(a[1],b[1],g),mul(a[2],b[2],g),g),mul(a[3],b[3],g),g))
end

"Return the bounded three-dimensional cross product as an `SVector`."
@inline function cross(a::StaticVector{3,CrossingFloat},
                       b::StaticVector{3,CrossingFloat})
    return SVector{3,CrossingFloat}(cross3(Tuple(a), Tuple(b)))
end

"Return the bounded three-dimensional dot product."
@inline dot(a::StaticVector{3,CrossingFloat}, b::StaticVector{3,CrossingFloat}) =
    dot3(Tuple(a), Tuple(b))

@inline cross(a::NTuple{3,CrossingFloat}, b::NTuple{3,CrossingFloat}) = cross3(a, b)
@inline dot(a::NTuple{3,CrossingFloat}, b::NTuple{3,CrossingFloat}) = dot3(a, b)


# An admitted vector of exact Float64 coordinates. The tokenized inner
# constructor keeps unchecked construction inside this module.
struct _Exact3Token end
const _EXACT3_TOKEN = _Exact3Token()
struct Exact3
    xyz::NTuple{3,Float64}
    Exact3(xyz::NTuple{3,Float64}, ::_Exact3Token) = new(xyz)
end

"""
    exact3(xyz::NTuple{3,Float64}) -> Union{Exact3,Nothing}

Admit three exact Float64 coordinates once for the specialized vector kernels.
"""
@inline function exact3(xyz::NTuple{3,Float64})
    admissible(xyz[1]) & admissible(xyz[2]) & admissible(xyz[3]) || return nothing
    return Exact3(xyz, _EXACT3_TOKEN)
end

@inline function _exact_difference_product(a::Float64, b::Float64,
                                           c::Float64, d::Float64)
    ph, pl = twoprod(a, b)
    qh, ql = twoprod(c, d)
    (ph == qh) & (pl == ql) && return raw(0.0, 0.0, 0.0)
    sh, se = twosum(ph, -qh)
    t = se + pl
    t2 = t - ql
    rh, rl = twosum(sh, t2)
    bound = E * ((abs(se) + abs(pl)) + (abs(t) + abs(ql)))
    return raw(rh, rl, pad(bound))
end

# Exact3 inputs make all input low limbs and radii statically zero. Each
# component pays one final admission after its two products are combined.
@inline function cross3(a::Exact3, b::Exact3)
    x = a.xyz
    y = b.xyz
    return (admit(_exact_difference_product(x[2], y[3], x[3], y[2])),
            admit(_exact_difference_product(x[3], y[1], x[1], y[3])),
            admit(_exact_difference_product(x[1], y[2], x[2], y[1])))
end

@inline function sum3(a::Exact3, b::Exact3)
    x = a.xyz
    y = b.xyz
    h1, l1 = twosum(x[1], y[1])
    h2, l2 = twosum(x[2], y[2])
    h3, l3 = twosum(x[3], y[3])
    return (make(h1, l1, 0.0), make(h2, l2, 0.0), make(h3, l3, 0.0))
end


end
