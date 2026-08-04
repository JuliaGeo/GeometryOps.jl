#=
# Sample

```@example sample
import GeoInterface as GI, GeometryOps as GO
p1 = GI.Polygon([[[-55965.680060140774, -31588.16072168928], [-55956.50771556479, -31478.09258677756], [-31577.548550575284, -6897.015828572996], [-15286.184961223798, -15386.952072224134], [-9074.387601621409, -27468.20712382156], [-8183.4538916097845, -31040.003969070774], [-27011.85123029944, -38229.02388009402], [-54954.72822634951, -32258.9734800704], [-55965.680060140774, -31588.16072168928]]])
points = GO.sample(p1, 100)
using CairoMakie
f, a, p = poly(p1)
scatter!(a, points)
f
```

=#

export sample, UniformSampling

struct UniformSampling
end

application_level(::UniformSampling) = TraitTarget(GI.MultiPolygonTrait(), GI.MultiLineStringTrait(), GI.MultiPointTrait(), GI.PolygonTrait(), GI.LineStringTrait())

function sample(geom, n::Int)
    return sample(UniformSampling(), geom, n)
end

function sample(alg, geom, n)
    return apply(x -> _sample(alg, GI.trait(x), x, n), application_level(alg), geom)
end

function _sample(alg::UniformSampling, ::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geom, n)
    (; X, Y) = GI.extent(geom)
    points = fill((0.0, 0.0), n)
    i = 1
    while i <= n
        x = rand() * (X[2] - X[1]) + X[1]
        y = rand() * (Y[2] - Y[1]) + Y[1]
        if contains(geom, (x, y))
            points[i] = (x, y)
            i += 1
        end
    end
    return points
end

function _sample(alg::UniformSampling, ::GI.LineStringTrait, geom, n)
    edges = to_edges(geom)
    edge_lengths = map(splat(distance), edges)
    # normalize the vector
    edge_probabilities = edge_lengths ./ sum(edge_lengths)
    edge_idxs = 1:length(edges)
    return map(1:n) do _
        edge_idx = sample(edge_idxs, edge_probabilities)
        x1, y1 = edges[edge_idx][1]
        x2, y2 = edges[edge_idx][2]
        distance = edge_lengths[edge_idx]
        t = rand() * distance
        (x1 + t * (x2 - x1), y1 + t * (y2 - y1))
    end
end