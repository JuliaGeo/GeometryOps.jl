# ValidatedFloats (internal module)

`ValidatedFloats` is a fixed two-limb `Float64` prototype for the narrow
arithmetic used by GeometryOps crossing emission. It is not a general interval
package or an arbitrary precision multifloat. The planar and spherical crossing emitters use the restricted domain.
The module lives inside GeometryOps and adds no dependencies beyond its existing
LinearAlgebra and StaticArrays dependencies.

Each `ValidatedFloat` represents an exact two-limb center `hi + lo` and a
nonnegative absolute radius `r`. While `isbounded(x)` is true, the represented
real is in `[hi + lo - r, hi + lo + r]`. Operations whose proof preconditions
cannot be established return an infinite radius, and `certify(Float64, x)` then
returns `nothing`.

The supported surface is exact construction from `Float64`, `+`, `-`, `*`,
`/`, `sqrt`, power-of-two `ldexp`, a strict Float64 rounding certificate, and
`certified_sign` for geometry decisions without inspecting the stored limbs.
The implementation rejects nonfinite arithmetic, unsafe exponent ranges,
inexact underflow, division by an interval containing zero, and malformed
two-limb centers.

## Restricted crossing domain

The nested `ValidatedFloats.CrossingFloats` module is a separate, faster domain
for the crossing emitter. `CrossingFloat` accepts only zero or finite leading
values with magnitude from `2^-400` through `2^400`; its division and square
root methods additionally require the normalized ranges documented in their
methods. Unsupported values produce an unbounded result whose certificate
fails.

Division accepts an admitted nonzero denominator when
`abs(y.lo) + y.rad ≤ 0.19abs(y.hi)`. Denominators outside the kernel's fast
range `1 ≤ abs(y.hi) ≤ 4` are normalized by scaling numerator and denominator
by the same exact power of two. If either scaled operand leaves the admitted
domain, division fails closed with an unbounded result. Restricted square root
requires `1 ≤ x.hi ≤ 16` and `abs(x.lo) + x.rad ≤ 0.19x.hi`.

The public vector interface uses `StaticArrays` and `LinearAlgebra`:

```julia
using LinearAlgebra, StaticArrays
using GeometryOps.ValidatedFloats.CrossingFloats

a = SVector(CrossingFloat.((0.8, 0.6, 0.0)))
b = SVector(CrossingFloat.((0.8, 0.0, 0.6)))
c = cross(a, b)                 # SVector{3,CrossingFloat}
d = dot(a, b)                   # CrossingFloat
s = a + b                       # SVector{3,CrossingFloat}
rounded = certify.(Ref(Float64), c)
```

Three-element tuples of `CrossingFloat` also support the fused `cross` and
`dot` methods; tuple `cross` returns a tuple. Other nonempty tuple lengths use
Julia's generic `dot` with validated scalar arithmetic.

The fused tuple kernels and exact-input `Exact3` experiment are private
implementation and audit helpers. The supported vector interface is
`LinearAlgebra.cross` and `dot` on static vectors and three-element tuples.
See [`docs/crossing-floats-proof.md`](docs/crossing-floats-proof.md) for the
local scalar/fused proof outline and test scope, and
[`docs/exact3-proof.md`](docs/exact3-proof.md) for the specialized kernel audit.

Both arithmetic domains remain internal and experimental. The module, tests,
and proof notes are kept together by topic to allow extraction into a separate
package later.

## General `ValidatedFloat` bound construction

The general `ValidatedFloat` addition and multiplication propagate input radii and outward-round every
nonnegative bound operation. TwoSum and TwoProduct preserve the exact leading
operation in their guarded normal ranges. The general multiplication bound can introduce a positive radius even for
exact Float64 inputs; the restricted `CrossingFloat` exact-input product path
preserves zero radius when its preconditions hold.

Division first computes a two-limb candidate `q`, then encloses its error using

```text
|x/y - q| = |x - q*y| / |y|.
```

The numerator is evaluated with validated arithmetic and the denominator uses
an outward lower bound for `|y|`. Square root similarly computes a candidate
`q` and uses

```text
|sqrt(x) - q| = |x - q^2| / (sqrt(x) + q)
                <= |x - q^2| / lower(q),  q > 0.
```

These identities remove expression-specific safety factors. The enclosure
arguments are documented and independently reviewed, supported by exact-rational
audits. A fallback handles failed certificates; it cannot repair an incorrect
successful certificate.

## Source layout

The module entry points include implementation files from [`src/`](src/):

- [`validated_float.jl`](src/validated_float.jl): general type and certificates.
- [`bounds.jl`](src/bounds.jl): error-free transforms and outward rounding.
- [`arithmetic.jl`](src/arithmetic.jl): general scalar operations.
- [`crossing_float.jl`](src/crossing_float.jl): restricted domain, scaling, and certificates.
- [`crossing_arithmetic.jl`](src/crossing_arithmetic.jl): restricted scalar operations.
- [`crossing_vectors.jl`](src/crossing_vectors.jl): fused `cross` and `dot` methods.
- [`exact_vectors.jl`](src/exact_vectors.jl): exact-input vector experiments.

## Running tests

From the GeometryOps repository root, run the internal module tests in the docs
workspace environment:

```bash
julia --project=docs -e 'push!(LOAD_PATH, abspath("test")); include("test/utils/ValidatedFloats/runtests.jl")'
```

For a persistent Julia daemon, the equivalent command is:

```bash
jld --project=docs eval 'push!(LOAD_PATH, abspath("test")); include("test/utils/ValidatedFloats/runtests.jl")'
```
