# Shared Natural Earth corpus for the overlay validation suites.
#
# Both `test/methods/clipping/overlayng/realdata_identities.jl` (oracle-free
# identities) and `test/methods/clipping/overlayng/s2_differential.jl`
# (differential against s2geography) sweep the same country polygons, so the
# loading, filtering and pairing rules live here once.
#
# Availability is the caller's business: `load_ne` throws if NaturalEarth /
# GeoJSON are missing, and every caller gates on that.

import GeometryOps as GO
import GeoInterface as GI
import LibGEOS as LG
import Extents

# Natural Earth 110 m **Sudan** is a documented bad input: its exterior ring is
# clockwise, carries a duplicate vertex, and is geodesically self-intersecting.
# s2geography rejects it too. It is excluded by name rather than left to poison
# a sweep — the engine contracts on valid input (design §2.2).
const NE_EXCLUDED = Set(["Sudan"])

"""
    load_ne(resolution) -> (names, geoms)

Natural Earth `admin_0_countries` at `resolution` (110, 50 or 10), as
`GO.tuples` polygons. Drops the documented bad inputs (`NE_EXCLUDED`) and any
ring GEOS calls invalid — the engine contracts on valid input.
"""
function load_ne(resolution)
    names = String[]
    geoms = Any[]
    fc = NaturalEarth.naturalearth("admin_0_countries", resolution)
    for f in fc
        g = GeoJSON.geometry(f)
        (g === nothing || GI.npoint(g) == 0) && continue
        nm = try string(f.NAME) catch; "?" end
        nm in NE_EXCLUDED && continue
        t = GI.trait(g)
        (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) || continue
        tg = GO.tuples(g)
        #-- NE has a handful of other self-intersecting rings besides Sudan
        try
            LG.isValid(GI.convert(LG, tg)) || continue
        catch
            continue
        end
        push!(names, nm); push!(geoms, tg)
    end
    return names, geoms
end

"""
    neighbour_pairs(names, geoms, limit) -> Vector{Tuple{Int,Int}}

Pairs whose extents intersect, in a deterministic order (sorted by name), up to
`limit`. Only pairs that actually share more than a point are kept, so a sweep
does not fill up with trivially-disjoint neighbours.
"""
function neighbour_pairs(names, geoms, limit)
    order = sortperm(names)
    exts = [GI.extent(g) for g in geoms]
    pairs = Tuple{Int, Int}[]
    for ii in eachindex(order), jj in (ii + 1):length(order)
        i, j = order[ii], order[jj]
        (exts[i] === nothing || exts[j] === nothing) && continue
        Extents.intersects(exts[i], exts[j]) || continue
        push!(pairs, (i, j))
        length(pairs) >= 4 * limit && break
    end
    kept = Tuple{Int, Int}[]
    for (i, j) in pairs
        length(kept) >= limit && break
        try
            LG.intersects(GI.convert(LG, geoms[i]), GI.convert(LG, geoms[j])) || continue
        catch
            continue
        end
        push!(kept, (i, j))
    end
    return kept
end

"""
    shifted(A, dx, dy)

`A` translated by `(dx, dy)` degrees. `A` against `shifted(A, ...)` makes every
border segment cross its own copy — the densest crossing workload this data has.
"""
shifted(A, dx, dy) = GO.apply(GI.PointTrait(), A) do p
    (GI.x(p) + dx, GI.y(p) + dy)
end

# Countries away from the antimeridian (whose overlay is a separate concern,
# handled by `antimeridian_split`), used for the shifted-self cases.
const NE_SHIFT_NAMES = ["Brazil", "France", "Egypt", "Australia", "Chile", "Norway",
                        "Indonesia", "India", "Kazakhstan", "Argentina"]

"""
    shifted_cases(names, geoms; dx = 0.5, dy = 0.25) -> Vector{Tuple{String,Any,Any}}

`(label, A, shifted(A))` for every `NE_SHIFT_NAMES` country present in `names`.
"""
function shifted_cases(names, geoms; dx = 0.5, dy = 0.25)
    out = Any[]
    for nm in NE_SHIFT_NAMES
        idx = findfirst(==(nm), names)
        idx === nothing && continue
        A = geoms[idx]
        push!(out, ("$nm shifted", A, shifted(A, dx, dy)))
    end
    return out
end
