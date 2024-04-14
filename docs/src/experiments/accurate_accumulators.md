# Accurate accumulation

Accurate arithmetic is a technique which allows you to calculate using more precision than the provided numeric type.

We will use the accurate sum routines from [AccurateArithmetic.jl](https://github.com/JuliaMath/AccurateArithmetic.jl) to show the difference!

```@example accurate
import GeometryOps as GO, GeoInterface as GI
using GeoJSON
using AccurateArithmetic
using NaturalEarth

all_adm0 = naturalearth("admin_0_countries", 10)
```
```@example accurate
GO.area(all_adm0)
```
```@example accurate
AccurateArithmetic.sum_oro(GO.area.(all_adm0.geometry))
```

```@example accurate
AccurateArithmetic.sum_kbn(GO.area.(all_adm0.geometry))
```

```@example accurate
GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |> collect .|> x -> [x]) .|> GO.signed_area |> sum
```

```@example accurate
GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |> collect .|> x -> [x]) .|> GO.signed_area |> sum_oro

```@example accurate
GI.Polygon.(GO.flatten(Union{GI.LineStringTrait, GI.LinearRingTrait}, all_adm0) |> collect .|> x -> [x]) .|> GO.signed_area |> sum_kbn
```