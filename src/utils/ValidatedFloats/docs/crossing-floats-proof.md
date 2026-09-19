# Restricted crossing arithmetic proof notes

`ValidatedFloats.CrossingFloats` is a specialized experimental arithmetic
domain for the planar and spherical crossing emitters. Its admitted nonzero
leading limbs are restricted to `[2^-400, 2^400]`. This keeps Float64
TwoProduct residuals representable throughout the intended determinant, cross,
dot, normalization, square-root, and division chains. Results outside the
admitted scalar domain become unbounded and cannot certify.

Every `CrossingFloat` represents `hi + lo ± rad`. TwoSum and TwoProduct supply
the exact leading operation where their guarded range assumptions hold. The
remaining rounded steps contribute explicit first-order bounds, padded by
`1 + 2^-44` and a `2^-1000` absolute floor. Division and square root are
restricted to the normalized denominator domains checked in their kernels.
Public division preserves the quotient while moving an arbitrary admitted,
well-separated denominator into that domain: it multiplies both operands by
the same power of two selected from the denominator's leading exponent.
`scale_pow2` admits the operation only when every nonzero limb and radius stays
representable and the scaled values stay in the arithmetic domain. Otherwise
the result is unbounded. Denominators already in `[1, 4]` take the original
kernel directly.

The private `Exact3` audit type admits an `NTuple{3,Float64}` once. Its inner
constructor requires a module-private token. The specialized private `cross3`
computes two exact Float64 products for each coordinate, combines
their expansions, attaches the same padded rounding bound as scalar addition,
then performs one final admission. Identical product expansions produce an
exact zero with zero radius. `sum3` uses exact TwoSum per coordinate and admits
the resulting expansion.

The local tests compare fused scalar and Exact3 vector kernels with exact
`Rational{BigInt}` coordinates, including coordinate-plane exact zeros. The
candidate also passed the external adversarial enclosure audit recorded during
the spike. These tests support the hand proof; they do not make it a
machine-checked proof.
