# Spherical Polygon Intersection - Root Cause Analysis

## Summary
- **Total failing cases**: 1008 (out of 1058 features)
- **Successful**: 50

## Error Type Breakdown (Updated 2026-01-10)
| Error Type | Count | Representative Index | Status |
|------------|-------|---------------------|--------|
| BoundsError | 779 | 1 | ROOT CAUSE IDENTIFIED |
| AssertionError | 225 | 2 | ROOT CAUSE IDENTIFIED |
| UndefVarError | ~~50~~ 0 | 7 | **RESOLVED** |
| TracingError | 4 | 317 | ROOT CAUSE IDENTIFIED |

---

## Root Cause Analysis

### 1. BoundsError (779 cases)
**Representative case**: Feature 1
**Status**: ROOT CAUSE IDENTIFIED

#### Full Stack Trace
```
BoundsError: attempt to access 5-element Vector{GeometryOps.PolyNode{Float64, GeometryOps.UnitSpherical.UnitSphericalPoint{Float64}}} at index [0]
Stacktrace:
  [1] throw_boundserror(A::Vector{...}, I::Tuple{Int64})
    @ Base ./essentials.jl:15
  [2] getindex
    @ ./essentials.jl:919 [inlined]
  [3] _classify_crossing!(alg::..., ::Type{Float64}, a_list::Vector{...}, b_list::Vector{...}; exact::...)
    @ GeometryOps src/methods/clipping/clipping_processor.jl:798
  [4] _classify_crossing!
    @ src/methods/clipping/clipping_processor.jl:711 [inlined]
  [5] _build_ab_list(alg::..., poly_a::..., poly_b::..., delay_cross_f::..., delay_bounce_f::...; exact::...)
    @ GeometryOps src/methods/clipping/clipping_processor.jl:176
  [6-12] ... intersection call chain ...
```

#### What Line/Function Causes the Error
The error occurs at **line 798** in `src/methods/clipping/clipping_processor.jl`:
```julia
start_pt = a_list[start_chain_idx]  # start_chain_idx is 0!
```

This is inside the `_classify_crossing!` function, in the block that handles closing an overlapping chain that started before the loop began (lines 785-807).

#### What Causes the Bad Index Value

The logic flow:
1. At line 719: `start_chain_idx` is initialized to `0`
2. At line 720: `unmatched_end_chain_edge` is initialized to `unknown`
3. The loop processes intersection points in order
4. When encountering an "end of overlapping chain" (line 744) AND `start_chain_edge == unknown` (line 746):
   - `unmatched_end_chain_edge` is set to the computed `b_side`
   - `unmatched_end_chain_idx` is set to `i`
5. The expectation is that later in the loop, we'll encounter the "start of overlapping chain" (line 773) which sets `start_chain_idx`
6. After the loop (line 785), if `unmatched_end_chain_edge != unknown`, the code assumes `start_chain_idx` was set

**THE BUG**: The code at line 785 only checks `unmatched_end_chain_edge != unknown`, but it doesn't verify that `start_chain_edge != unknown` (or equivalently, that `start_chain_idx != 0`). When we only encounter the "end" of an overlapping chain but never the "start", `start_chain_idx` remains 0.

#### Geometry That Triggers This
The test case has two polygons with a **shared edge segment**:
- **Poly1**: Rectangle from (70,-81) to (71,-80)
  - Has edge: (70,-81) -> (70,-80)
- **Poly2**: Rectangle from (69,-80.47813) to (70,-79.52949)
  - Has edge: (70,-79.52949) -> (70,-80.47813)

These edges **overlap** along x=70 from y=-80 to y=-80.47813. This creates an "overlapping chain" in the Foster-Hormann clipping algorithm.

The problem occurs because:
1. The algorithm starts processing from the middle of an overlapping chain
2. It encounters the "end" of the chain first (setting `unmatched_end_chain_edge`)
3. But due to the geometry, it never encounters a corresponding "start" of the chain
4. The post-loop cleanup at line 785 runs (because `unmatched_end_chain_edge != unknown`)
5. It tries to access `a_list[0]` because `start_chain_idx` was never updated from its initial value of 0

#### Root Cause Hypothesis

**The root cause is a missing guard condition at line 785.**

The condition `if unmatched_end_chain_edge != unknown` is insufficient. It should also verify that `start_chain_edge != unknown` (or `start_chain_idx != 0`) before attempting to access `a_list[start_chain_idx]`.

Possible fixes:
1. Add a guard: `if unmatched_end_chain_edge != unknown && start_chain_edge != unknown`
2. Or handle the degenerate case where we have an unmatched end but no start (this may indicate the entire ring is an overlapping chain, requiring special handling)

**Why this affects Spherical more than Planar**: Testing confirms that `Planar()` succeeds on these exact same polygons (returning 0 result polygons, which is correct since they only share an edge with no overlapping area). The spherical geometry calculations in `_get_sides()` likely produce different orientation results that cause the "end of overlapping chain" condition to be triggered when it shouldn't be, or fail to trigger the "start of overlapping chain" condition when it should be.

The key difference is in how `_get_sides()` computes which side of a line segment a point is on - this uses spherical orientation predicates via `_is_collinear()` and `_side_of()` which may give different results than the planar versions for near-collinear points near the poles (these polygons are at latitude -80 to -81, close to the South Pole).

**Verification**: Running `GO.intersection(GO.Planar(), poly1, poly2)` succeeds, while `GO.intersection(GO.Spherical(), poly1, poly2)` fails with BoundsError.

#### Deeper Analysis: Why Spherical Triggers This

The key functions involved:
1. `_get_sides()` (line 851) - determines which side of the a-polygon edges the b-polygon neighbors are on
2. `_get_side()` (line 901 for Spherical) - uses `spherical_orient()` to compute orientation
3. `spherical_orient()` in `src/utils/UnitSpherical/predicates.jl` - uses cross product and dot product

The spherical predicates use a tolerance check:
```julia
tol = eps(Float64) * 16  # ~3.5e-15
if abs(dot_product) < tol
    return 0  # collinear
end
```

Near the poles (latitude -80 to -81), the coordinate transformations can amplify numerical errors. The `UnitSphereFromGeographic()` transformation converts lat/lon to 3D unit sphere coordinates:
- At latitude -80, points are very close to the south pole (z ≈ -0.98)
- Small longitude differences become very small xy differences
- This can cause orientation predicates to incorrectly return 0 (collinear) or flip signs

When these predicates return unexpected values, the `_classify_crossing!` logic misclassifies the relationship between polygon edges:
1. It may detect an "end of overlapping chain" when edges aren't actually overlapping
2. It may fail to detect the "start of overlapping chain"
3. This leaves `start_chain_idx` at its initial value of 0

#### Related File Locations
- Error location: `src/methods/clipping/clipping_processor.jl:798`
- Chain classification logic: `src/methods/clipping/clipping_processor.jl:711-807`
- Spherical side computation: `src/methods/clipping/clipping_processor.jl:901-919`
- Spherical orientation predicate: `src/utils/UnitSpherical/predicates.jl:30-43`

---

### 2. AssertionError (225 cases)
**Representative case**: Feature 2
**Status**: ROOT CAUSE IDENTIFIED
**Root cause**: Float32 precision too low for unit length tolerance check

#### Full Stack Trace
```
AssertionError: Input vector 'a' must be unit length
Stacktrace:
  [1] robust_cross_product(a::UnitSphericalPoint{Float32}, b::UnitSphericalPoint{Float32})
    @ GeometryOps.UnitSpherical.RobustCrossProduct src/utils/UnitSpherical/robustcrossproduct/RobustCrossProduct.jl:99
  [2] spherical_orient(a::UnitSphericalPoint{Float32}, b::UnitSphericalPoint{Float32}, c::UnitSphericalPoint{Float64})
    @ GeometryOps.UnitSpherical src/utils/UnitSpherical/predicates.jl:33
  [3] point_on_spherical_arc(p::UnitSphericalPoint{Float64}, a::UnitSphericalPoint{Float32}, b::UnitSphericalPoint{Float32})
    @ GeometryOps.UnitSpherical src/utils/UnitSpherical/predicates.jl:75
  [4] _point_filled_curve_orientation(::Spherical{Float64}, point::UnitSphericalPoint{Float64}, curve::GeoJSON.LineString{2, Float32}; ...)
    @ GeometryOps src/methods/geom_relations/geom_geom_processors.jl:596
  [5] _pt_off_edge_status @ src/methods/clipping/clipping_processor.jl:940
  [6] _flag_ent_exit! @ src/methods/clipping/clipping_processor.jl:969
  [7] _build_ab_list @ src/methods/clipping/clipping_processor.jl:180
  [8] _intersection @ src/methods/clipping/intersection.jl:84
  [9] intersection @ src/methods/clipping/intersection.jl:42
```

#### What Assertion is Failing
The assertion at `RobustCrossProduct.jl:99`:
```julia
@boundscheck @assert isUnitLength(a) "Input vector 'a' must be unit length"
```

The `isUnitLength` function in `robustcrossproduct/utils.jl:28-30`:
```julia
function isUnitLength(v::AbstractVector)
    return isapprox(sum(abs2, v), 1.0, rtol=1e-14)
end
```

#### What Condition Causes the Failure
The tolerance `rtol=1e-14` is appropriate for Float64 but far too strict for Float32:
- **Float64 machine epsilon**: ~2.2e-16 --> tolerance of 1e-14 is ~100x epsilon (reasonable)
- **Float32 machine epsilon**: ~1.2e-7 --> tolerance of 1e-14 is ~10,000,000x smaller than epsilon (impossible to satisfy)

When geographic coordinates (longitude, latitude) are converted to UnitSphericalPoint via `sincosd()`:
- Float32 coordinates produce Float32 UnitSphericalPoints
- The computed norm squared deviates from 1.0 by ~1e-7 to 1e-8 (on the order of Float32 epsilon)
- Example: norm squared = 0.99999994 (deviation = 6e-8) or norm squared = 1.0000001 (deviation = 1.2e-7)

**Data flow causing the issue**:
1. GeoJSON file stores coordinates as Float32 (`GeoJSON.Point{2, Float32}`)
2. `UnitSphereFromGeographic()` converts to `UnitSphericalPoint{Float32}` (preserves type)
3. `point_on_spherical_arc` calls `spherical_orient` with Float32 points
4. `spherical_orient` calls `robust_cross_product`
5. `robust_cross_product` asserts `isUnitLength(a)` which fails for Float32

#### Root Cause Hypothesis
**Primary cause**: The `isUnitLength` function uses a fixed tolerance of `1e-14` that does not account for the precision of the input type. Float32 inputs cannot satisfy this tolerance due to fundamental precision limitations.

**Contributing factors**:
1. GeoJSON.jl stores coordinates as Float32 by default
2. Type preservation in `UnitSphericalPoint` constructors propagates Float32 through the pipeline
3. No type promotion to Float64 occurs before the assertion check

#### Potential Fixes
1. **Make `isUnitLength` tolerance type-aware**: Use `rtol = 16 * eps(eltype(v))` instead of fixed `1e-14`
2. **Promote Float32 to Float64 in coordinate conversion**: Ensure `UnitSphereFromGeographic` always returns Float64
3. **Normalize after conversion**: Add explicit normalization step after geographic-to-spherical conversion
4. **Both 1 and 2**: Most robust solution - adapt tolerance AND ensure Float64 precision in spherical calculations

#### Related File Locations
- Assertion that fails: `src/utils/UnitSpherical/robustcrossproduct/RobustCrossProduct.jl:99`
- isUnitLength function: `src/utils/UnitSpherical/robustcrossproduct/utils.jl:28-30`
- Coordinate transform: `src/utils/UnitSpherical/coordinate_transforms.jl:38-56`
- UnitSphericalPoint constructors: `src/utils/UnitSpherical/point.jl:41-76`

---

### 3. UndefVarError (Previously 50 cases, NOW 0)
**Representative case**: Feature 7
**Status**: RESOLVED - No longer occurring

#### Investigation Result
Feature 7 now succeeds, returning an empty polygon array:
```
Feature 7:
  Poly1 points: 5
  Poly2 points: 5
  SUCCESS: 0 polygons
```

Re-running the full test suite confirms that **UndefVarError no longer occurs for any feature**.

Updated error counts after re-investigation:
| Error Type | Old Count | Current Count |
|------------|-----------|---------------|
| BoundsError | 779 | 779 |
| AssertionError | 225 | 225 |
| UndefVarError | 50 | **0** |
| TracingError | 4 | 4 |
| Success | 0 | **50** |

#### Likely Fix
The UndefVarError was likely fixed in one of these recent commits:
- `bd37ae516 Fix GI.getcoord on UnitSphericalPoint`
- `faa0f74bc Fix slerp returning inf/nan values for equal points`

These commits addressed edge cases in spherical point handling that would have caused undefined variable errors in downstream calculations.

---

### 4. TracingError (4 cases)
**Representative case**: Feature 317
**Status**: ROOT CAUSE IDENTIFIED
**Root cause**: Duplicate vertex at pole creates inconsistent neighbor relationships in polygon tracing

#### Full Stack Trace
```
TracingError: Clipping tracing hit every point - clipping error.
Please open an issue with the polygons contained in this error object.

Polygon A:
Vector{Tuple{Float32, Float32}}[[(159.0, 89.0), (180.0, 90.0), (180.0, 90.0), (160.0, 89.0), (159.0, 89.0)]]

Polygon B:
Vector{Tuple{Float32, Float32}}[[(-20.0, 89.52514), (160.0, 89.52514), (106.975914, 89.21064), (33.02409, 89.21064), (-20.0, 89.52514)]]

Stacktrace:
  [1] _trace_polynodes(alg::GeometryOps.FosterHormannClipping{GeometryOpsCore.Spherical{Float64}, GeometryOps.NestedLoop}, ::Type{Float64}, a_list::Vector{GeometryOps.PolyNode{Float64, GeometryOps.UnitSpherical.UnitSphericalPoint{Float64}}}, b_list::Vector{GeometryOps.PolyNode{Float64, GeometryOps.UnitSpherical.UnitSphericalPoint{Float64}}}, a_idx_list::Vector{Int64}, f_step::typeof(GeometryOps._inter_step), poly_a::GeoJSON.Polygon{2, Float32}, poly_b::GeoJSON.Polygon{2, Float32})
    @ GeometryOps ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/clipping_processor.jl:1093
  [2] _intersection(alg::GeometryOps.FosterHormannClipping{GeometryOpsCore.Spherical{Float64}, GeometryOps.NestedLoop}, ::GeometryOpsCore.TraitTarget{GeoInterface.PolygonTrait}, ::Type{Float64}, ::GeoInterface.PolygonTrait, poly_a::GeoJSON.Polygon{2, Float32}, ::GeoInterface.PolygonTrait, poly_b::GeoJSON.Polygon{2, Float32}; exact::GeometryOpsCore.True, kwargs::@Kwargs{})
    @ GeometryOps ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:85
  [3] _intersection
    @ ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:74 [inlined]
  [4] #intersection#150
    @ ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:45 [inlined]
  [5] intersection (repeats 2 times)
    @ ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:42 [inlined]
  [6] #intersection#152
    @ ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:59 [inlined]
  [7] intersection
    @ ~/.julia/dev/geo/GeometryOps.jl/src/methods/clipping/intersection.jl:58 [inlined]
```

#### What Causes the Tracing to Fail

The `_trace_polynodes` function at line 1093 throws a `TracingError` when `visited_pts >= total_pts`. This happens because the tracing algorithm gets stuck in an infinite loop, cycling between the same points without ever closing the polygon.

**The infinite loop pattern:**

The tracing state shows the algorithm repeatedly cycles through:
1. `a_list[4]` → neighbor is `b_list[4]` → switch to b_list
2. `b_list[4]` → step forward to `b_list[3]` then `b_list[4]` → neighbor is `a_list[4]` → switch to a_list
3. `a_list[4]` → step backward to `a_list[3]` → neighbor is `b_list[2]` → switch to b_list
4. `b_list[2]` → step forward to `b_list[3]` then `b_list[4]` → back to step 2

This happens because:
- `b_list[3]` and `b_list[4]` have **essentially the same point** (distance ~8.67e-19, i.e., floating point error)
- But `b_list[3]` is NOT marked as an intersection (`inter=false`)
- While `b_list[4]` IS marked as an intersection (`inter=true, crossing=true`)

The tracing algorithm expects to close the polygon by returning to the start point or its b_list neighbor, but the inconsistent point/intersection structure prevents this.

#### What Condition in the Polygon Geometry Leads to This

**The input geometry:**

Polygon A (poly1):
```
Point 1: (159.0, 89.0)
Point 2: (180.0, 90.0)  ← NORTH POLE
Point 3: (180.0, 90.0)  ← NORTH POLE (DUPLICATE!)
Point 4: (160.0, 89.0)
Point 5: (159.0, 89.0)  ← closing point
```

Polygon B (poly2):
```
Point 1: (-20.0, 89.52514)
Point 2: (160.0, 89.52514)
Point 3: (106.975914, 89.21064)
Point 4: (33.02409, 89.21064)
Point 5: (-20.0, 89.52514)  ← closing point
```

**Key observations:**

1. **Polygon A has a duplicate vertex at the North Pole** (180°, 90°). Points 2 and 3 are identical.

2. In UnitSpherical coordinates (3D unit sphere), both points become exactly `(0, 0, 1)` - the same point regardless of longitude because all longitudes converge at the pole.

3. This creates a **degenerate zero-length edge** from point 2 to point 3 in polygon A.

4. When the intersection algorithm processes this:
   - Edge 1-2 of Poly A creates an intersection point (at the pole)
   - Edge 2-3 of Poly A is degenerate (zero length)
   - Edge 3-4 of Poly A creates another intersection point very close to the pole

5. The intersection detection creates entries in `b_list` where:
   - `b_list[3]` is a regular vertex of poly_b (point 2→3 edge start)
   - `b_list[4]` is an intersection point
   - These two points are essentially at the same location but have different intersection status

#### The a_list and b_list State

```
a_list (5 nodes):
  [1] non-intersection vertex
  [2] crossing intersection, neighbor=5 (b_list[5]), EXIT
  [3] crossing intersection, neighbor=2 (b_list[2]), ENTER  ← AT NORTH POLE
  [4] crossing intersection, neighbor=4 (b_list[4]), EXIT
  [5] non-intersection vertex

b_list (7 nodes):
  [1] non-intersection vertex
  [2] crossing intersection, neighbor=3 (a_list[3]), ENTER  ← AT NORTH POLE
  [3] non-intersection vertex  ← nearly same point as b[4]!
  [4] crossing intersection, neighbor=4 (a_list[4]), EXIT   ← nearly same point as b[3]!
  [5] crossing intersection, neighbor=2 (a_list[2]), ENTER
  [6] non-intersection vertex
  [7] non-intersection vertex
```

The distance between `b_list[3].point` and `b_list[4].point` is ~8.67e-19 (essentially zero).

#### Root Cause Hypothesis

**Primary cause:** The input polygon has a **duplicate vertex at the geographic pole**, which creates a zero-length edge. This is a degenerate geometry that the algorithm does not handle correctly.

**Secondary cause:** The intersection detection algorithm does not consolidate nearly-identical points. When processing different edges, it creates separate entries for what are essentially the same point, leading to inconsistent neighbor relationships.

**Why tracing fails:**
1. The tracing algorithm follows edges and switches lists at crossing points
2. The neighbor relationship links `a_list[3]` → `b_list[2]` and `a_list[4]` → `b_list[4]`
3. But `b_list[3]` (between b_list[2] and b_list[4]) is not marked as an intersection
4. So when tracing from b_list[2] forward, it steps through b_list[3] (not a crossing) to b_list[4]
5. From b_list[4], it switches to a_list[4], then steps to a_list[3]
6. From a_list[3], it switches to b_list[2], repeating the cycle
7. The polygon never closes because the topology is inconsistent

**Contributing factors:**
1. Geographic coordinates converge at poles - all longitudes map to the same 3D point at lat=90°
2. No input validation to detect/remove duplicate vertices
3. No consolidation of nearly-identical intersection points
4. The Foster-Hormann algorithm assumes well-formed, non-degenerate input

#### Potential Fixes

1. **Input validation**: Detect and remove duplicate/nearly-identical vertices before processing
2. **Zero-length edge handling**: Skip edges where start == end in `_build_a_list`
3. **Point consolidation**: Merge nearly-identical points in b_list during `_build_b_list`
4. **Pole-aware processing**: Special handling for polygons that include/touch geographic poles
5. **Degeneracy detection**: Check for degenerate cases before tracing and handle appropriately

#### Related File Locations
- Error thrown: `src/methods/clipping/clipping_processor.jl:1093`
- Tracing algorithm: `src/methods/clipping/clipping_processor.jl:1055-1126`
- a_list construction: `src/methods/clipping/clipping_processor.jl:489-635`
- b_list construction: `src/methods/clipping/clipping_processor.jl:647-698`
- Spherical arc intersection: `src/utils/UnitSpherical/arc_intersection.jl:65-162`

---

## Summary and Recommended Fixes

### Priority Order

| Priority | Error Type | Count | Fix | Complexity |
|----------|------------|-------|-----|------------|
| 1 | BoundsError | 779 | Add `start_chain_edge != unknown` guard at line 785 | Low |
| 2 | AssertionError | 225 | Make `isUnitLength` tolerance type-aware: `rtol = 16 * eps(eltype(v))` | Low |
| 3 | TracingError | 4 | Pre-process to remove duplicate consecutive points | Medium |
| - | UndefVarError | 0 | Already fixed in recent commits | N/A |

### Key Code Locations

1. **BoundsError fix location**: `src/methods/clipping/clipping_processor.jl:785`
   - Change: `if unmatched_end_chain_edge != unknown`
   - To: `if unmatched_end_chain_edge != unknown && start_chain_edge != unknown`

2. **AssertionError fix location**: `src/utils/UnitSpherical/robustcrossproduct/utils.jl:29`
   - Change: `isapprox(sum(abs2, v), 1.0, rtol=1e-14)`
   - To: `isapprox(sum(abs2, v), 1.0, rtol = 16 * eps(eltype(v)))`

3. **TracingError fix location**: Input validation in `_build_a_list` or pre-processing step
   - Add: Remove consecutive duplicate points before processing

### Overall Assessment

The spherical polygon intersection has three main issues:

1. **Logic bug** (BoundsError): A missing guard condition allows indexing with 0 when chain processing is incomplete. This is a straightforward fix.

2. **Type precision mismatch** (AssertionError): The unit length check uses Float64-appropriate tolerance for Float32 inputs. This is also a straightforward fix.

3. **Degenerate geometry handling** (TracingError): Duplicate vertices at poles create inconsistent topology. This requires input validation or geometry pre-processing.

All three root causes have been identified with specific fix locations.
