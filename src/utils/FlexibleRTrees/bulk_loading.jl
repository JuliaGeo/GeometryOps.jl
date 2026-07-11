# # Bulk loading

# ## Leaf ordering

"""
    loadorder(algorithm, extents::Vector{<:Extents.Extent}, nodecapacity)::Vector{Int}

The permutation in which `algorithm` packs `extents` into leaves.  Implement
this for a new `BulkLoadAlgorithm` subtype to plug in another ordering.
"""
loadorder(::Unsorted, extents, nodecapacity) = collect(1:length(extents))
loadorder(::HPR, extents, nodecapacity) = sortperm(_hilbert_keys(extents))
function loadorder(::STR, extents::Vector{E}, nodecapacity) where E
    centers = [_center(e) for e in extents]
    perm = collect(1:length(extents))
    _str_tile!(perm, centers, 1, length(extents), 1, _ndims(E), nodecapacity)
    return perm
end

#=
One recursion level of N-dimensional sort-tile-recursive: sort the range by
the current dimension's center, cut it into `S ≈ P^(1/remaining)` slabs of
whole leaf pages, and tile each slab along the remaining dimensions.  After
the last dimension the consecutive `nodecapacity`-runs are the leaf tiles.
=#
function _str_tile!(perm, centers, lo, hi, dim, ndims, nodecapacity)
    len = hi - lo + 1
    len <= nodecapacity && return # a single leaf: internal order doesn't matter
    sort!(view(perm, lo:hi); by = i -> @inbounds(centers[i][dim]))
    dim == ndims && return
    P = cld(len, nodecapacity)                        # leaf pages in this range
    S = ceil(Int, P^(1 / (ndims - dim + 1)))          # slabs along this dimension
    slab = cld(P, S) * nodecapacity                   # items per slab (whole pages)
    i = lo
    while i <= hi
        _str_tile!(perm, centers, i, min(i + slab - 1, hi), dim + 1, ndims, nodecapacity)
        i += slab
    end
    return
end

# ## Packing

# Union consecutive runs of `nodecapacity` extents, bottom-up, until a level
# fits in one (implicit) root node.  Returns levels coarsest-first.
function _pack_levels(leaves::Vector{E}, nodecapacity::Int) where E
    levels = [leaves]
    current = leaves
    while length(current) > nodecapacity
        nparents = cld(length(current), nodecapacity)
        parents = Vector{E}(undef, nparents)
        for p in 1:nparents
            lo = (p - 1) * nodecapacity + 1
            hi = min(p * nodecapacity, length(current))
            acc = current[lo]
            for j in (lo + 1):hi
                acc = Extents.union(acc, @inbounds current[j])
            end
            parents[p] = acc
        end
        push!(levels, parents)
        current = parents
    end
    return reverse!(levels)
end

# ## Extent helpers

_ndims(::Type{Extents.Extent{K, V}}) where {K, V} = length(K)
_center(ext::Extents.Extent) = map(b -> (b[1] + b[2]) / 2, values(ext))
