# # Hilbert curve encoding
#
#=
The HPR ("Hilbert-packed R-tree") bulk loader sorts items by the position of
their center on a Hilbert space-filling curve.  The Hilbert curve visits every
cell of a `2^bits × … × 2^bits` grid exactly once, and consecutive cells along
the curve are always spatially adjacent — so items that are consecutive in
sorted order are close in space, in every dimension, at every scale.  That
fractal locality is what makes simple "pack consecutive runs" tree
construction produce tight nodes on all levels.

This is Skilling's transpose-based algorithm ("Programming the Hilbert
curve", AIP Conf. Proc. 707, 2004), which works in any number of dimensions —
unlike the lookup-table encoders (e.g. JTS's 2D-only `HilbertCode`).  The
"transpose" is the Hilbert index of the point stored bit-interleaved across
the N input coordinates; we finish by de-interleaving it into a single
integer sort key.
=#

"""
    hilbert_key(coords::NTuple{N, UInt32}, bits::Int) -> UInt64

The Hilbert-curve index of a point on the `N`-dimensional `2^bits` grid, as a
sortable integer.  Requires `N * bits <= 64`; each coordinate must be
`< 2^bits`.
"""
function hilbert_key(coords::NTuple{N, UInt32}, bits::Int) where N
    X = MVector{N, UInt32}(coords)
    M = one(UInt32) << (bits - 1)
    # Inverse undo
    Q = M
    while Q > one(UInt32)
        P = Q - one(UInt32)
        for i in 1:N
            if !iszero(X[i] & Q)
                X[1] ⊻= P # invert
            else # exchange
                t = (X[1] ⊻ X[i]) & P
                X[1] ⊻= t
                X[i] ⊻= t
            end
        end
        Q >>= 1
    end
    # Gray encode
    for i in 2:N
        X[i] ⊻= X[i - 1]
    end
    t = zero(UInt32)
    Q = M
    while Q > one(UInt32)
        !iszero(X[N] & Q) && (t ⊻= Q - one(UInt32))
        Q >>= 1
    end
    for i in 1:N
        X[i] ⊻= t
    end
    # Interleave the transpose into one key, most significant bits first.
    key = zero(UInt64)
    for b in (bits - 1):-1:0
        for i in 1:N
            key = (key << 1) | UInt64((X[i] >> b) & one(UInt32))
        end
    end
    return key
end

# Bits per dimension for the quantization grid: 16 where possible (the JTS
# HPRtree uses 12), fewer in high dimensions so the key fits in 64 bits.
hilbert_bits(N::Int) = min(16, 63 ÷ N)

#=
Quantize the centers of a set of extents onto the Hilbert grid spanned by
their total extent, and return each center's Hilbert key.  Degenerate
dimensions (zero span) collapse to grid coordinate 0.
=#
function _hilbert_keys(extents::Vector{E}) where E <: Extents.Extent
    N = _ndims(E)
    bits = hilbert_bits(N)
    total = reduce(Extents.union, extents)
    los = map(first, values(total))
    spans = map(b -> Float64(b[2] - b[1]), values(total))
    scale = Float64((1 << bits) - 1)
    keys = Vector{UInt64}(undef, length(extents))
    for (j, e) in enumerate(extents)
        c = _center(e)
        q = ntuple(Val(N)) do i
            span = spans[i]
            frac = iszero(span) ? 0.0 : (Float64(c[i]) - Float64(los[i])) / span
            UInt32(clamp(round(Int, frac * scale), 0, Int(scale)))
        end
        keys[j] = hilbert_key(q, bits)
    end
    return keys
end
