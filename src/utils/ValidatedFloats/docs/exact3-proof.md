# Exact3 specialized-kernel audit

`exact3` admits each exact Float64 coordinate once. Its sign-masked bit-range
test accepts signed zero or a finite magnitude in `[2^-400, 2^400]` and rejects
smaller nonzero values, larger values, infinities, and NaNs. The tokenized
constructor prevents unchecked public construction.

For `_exact_difference_product(a,b,c,d)`, TwoProduct gives exact expansions

`a*b = ph + pl` and `c*d = qh + ql`.

The admission range keeps both products finite and their nonzero FMA residuals
normal. Therefore `(ph == qh) && (pl == ql)` proves equality of the two exact
products, including signed permutations, and returning exact zero is sound.

Otherwise TwoSum gives `ph - qh = sh + se` exactly. The implementation then
forms

`t = RN(se + pl)` and `t2 = RN(t - ql)`

before an exact TwoSum of `sh + t2`. The only unrepresented errors are those
two ordinary additions. Binary64 unit roundoff bounds them by

`u * (|se| + |pl|) + u * (|t| + |ql|)`.

This is the implemented analytic bound with `E = u`. Its final FLOOR/PAD step
covers evaluation rounding and underflow exactly as in the audited scalar
arithmetic. Final `admit` either restores the public CrossingFloat invariant or
fails closed.

`cross3(::Exact3, ::Exact3)` applies this construction independently to each
coordinate. No unchecked product result becomes an operand of another product.

`sum3` applies TwoSum directly to each pair of exact Float64 coordinates, so
each returned two-limb center is exact with zero radius. A sum outside the
public leading-limb domain fails closed through `make`.

The module defines specialized Exact3 `cross3` and `sum3`; it does not define a
specialized `dot3(::Exact3, ::Exact3)`. The existing bounded CrossingFloat dot
kernel is unchanged. This is only an API observation if an exact-input dot
specialization was intended.
