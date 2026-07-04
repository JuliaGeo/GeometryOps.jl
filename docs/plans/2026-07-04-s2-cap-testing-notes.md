# S2 testing methodology notes — caps, bounds, exact predicates (2026-07-04)

Extracted from the s2geometry reference checkout
(`/Users/anshul/temp/GO_jts/s2geometry/src/s2/`, cites are `file:line`) to
guide testing of GeometryOps' unit-spherical cap/arc extents and robust cap
predicates.  Companion to
`2026-07-04-unit-spherical-indexing-foundations.md`.

## Cap predicate testing patterns

- Named fixtures (`s2cap_test.cc:55-170`): `empty`, `full`, axis singleton
  caps (`FromPoint`), `tiny` (radius 1e-10), `hemi`, `concave` (150°).
  Global `kEps = 1e-15` ("about 9× the double-precision roundoff relative
  error", `:52-53`).
- Tangent-radius accuracy (`:118-122`): perturb a tiny cap's center along a
  unit tangent by `0.99 r` (must contain) vs `1.01 r` (must not).
- **Bracketing/"sandwich" tolerance pattern** (`:133-148`): derive
  `max_error` per cap from the chord-angle constructor error methods, build
  `cap.PlusError(±max_error)`, and assert a boundary point is inside the
  expanded cap and outside the shrunk one.  Point-containment error grows
  with cap angle, so tolerance is per-cap, not a constant.
- Exact tangency (`:161-165`): `hemi.Contains(cap(π/4 − ε))` true,
  `π/4 + ε` false.
- Metamorphic laws: Union/Contains consistency (`:343-389`), `Expanded`
  monotonicity (`:314-323`), complement lattice (`:110-113`), and
  `ApproxEquals` with center-angle + radius² tolerance (`s2cap.cc:312-322`)
  as the canonical approximate cap equality.
- Cap-vs-cell tests use a **dot-product oracle** as cheap analytic truth
  (`:272-299`): the predicate under test is compared against a hand-derived
  dot-product threshold.

## Bound conservativeness + tightness (rect bounder)

- The tolerance is a centrally derived error budget:
  `MaxErrorForTests() = (10 ε lat, 1 ε lng)` with a per-term breakdown in
  comments (`s2latlng_rect_bounder.cc:345-357`).  Internal constants
  (`n_norm < 8.618 ε`, `m_error`, per-stage `+3ε`/`+2ε`/`+9ε` pads) each
  carry a justifying comment (`:67-74`, `:117-125`, `:143-147`, `:210`,
  `:339`).
- **Two-sided bound testing** (`NearlyIdenticalOrAntipodalPoints`,
  `s2latlng_rect_bounder_test.cc:179-244`, 10k iters):
  `bound ⊇ endpoints-box` (conservative) AND
  `endpoints-box.Expanded(kRectError) ⊇ bound` (tight).
- Adversarial generator `PerturbATowardsB` (`:136-161`): 10% identical
  point, 20% exactly-proportional (non-unit) point, 20% distance-squared
  underflow, else a **log-uniform ULP sweep** (`ε × LogUniform(1e-5, 10)`),
  dispatched across near-pole/near-equator regimes.
- Tightness anchored analytically: interior-extremum edges assert the
  result sits at the **middle of the error band**
  (`EXPECT_DOUBLE_EQ(truth + 0.5 kRectError, hi)`, `:83-84`); random edges
  are built **through a known extremum point** so truth is by construction
  (`MaxLatitudeRandom`, `:100-134`).
- Error-constant tightness probed by threshold-straddling test pairs
  (`:270-271`); regressions pinned with a `Sign`-predicate oracle
  (`AccuracyBug`, `:335-356`).

## Exact/filtered predicate testing

- 4-tier ladder (double → long double → exact → symbolic).  The **nesting
  invariant** (`:604-636`): any faster tier that returns non-zero must
  agree with every slower tier; returning 0 (uncertain) is the only
  allowed "failure".  Tests also assert *which* tier decides
  (`expected_prec`, `:626-629`).
- Coverage tests hand-craft inputs per tier via perturbation magnitude:
  `1e-15` → double, `1e-8` → long double, `1e-100` → exact, exact zeros →
  symbolic (`:642-659`); `SignDotProd` examples `(ε,1,0)`, `(1e-45,1,0)`
  (`:1177-1201`).
- Consistency stress (5000 iters, `:722-774`): `ChoosePoint` crushes each
  coordinate onto planes/axes with probability 1/3 by `LogUniform(1e-50,1)`
  (`:561-569`); ties are constructed **exactly on the decision boundary**
  via symmetric `GetPointOnLine` offsets (`:742-749`); a precision
  histogram gives tier-coverage visibility (`:533-556`).
- Degeneracy menu for orientation tests (`AddDegeneracy`, `:268-321`):
  points on the circle, 1-ulp and 1e-15 perturbations, `(1±ε)` rescales,
  coordinate-plane intersections, exactly-collinear tangent pairs,
  negations; hardest circles are the coordinate planes (`:372-391`).
- Filter give-up-rate is *quantified*, not just tolerated
  (`StableSignTest.FailureRate`, `:395-436`).

## Random fixture design

- Unit points: uniform-in-cube then normalize (`s2random.cc:47-55`) — not
  area-uniform, used pervasively anyway.
- Caps: **log-uniform in area** between bounds (`s2random.cc:92-105`,
  rationale `s2random.h:56-62`) — even coverage across scales, where
  error behavior differs.
- Area-uniform sampling inside a cap: uniform height `h`, uniform angle,
  circle radius `√(h(2−h))` through the cap's frame (`:107-128`).
- Determinism: per-test tagged seed hashed with one global seed flag
  (`s2testing.cc:157-163`); iteration counts are flags.

## Adopted in GeometryOps (test/utils/cap_extents.jl)

Two-sided per-axis tightness for `Extents.extent(cap)`; adversarial
axis-crushed points (`adversarial_point`); log-uniform radii; ±ulp
tangency sweeps solved in 512-bit precision against an independent
angle-space ground truth; exact-tie cap–point discrimination; consistency
laws incl. the not-quite-full sharp edge (`k = −1` with `‖c‖ > 1`);
explicit arc-midpoint (max sagitta) and degenerate-endpoint arc tests.
Skipped: S2Cell machinery, the long-double middle tier, symbolic
perturbation (our predicates resolve true ties exactly instead).
