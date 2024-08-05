#=
# Lazy wrappers

These wrappers lazily apply some fixes like closing rings.
=#

abstract type AbstractLazyWrapper{GeomType} end

struct LazyClosedRing{GeomType} <: AbstractLazyWrapper{GeomType}
    ring::GeomType
    function LazyClosedRing(ring)
        LazyClosedRing(GI.trait(ring), ring)
    end
    function LazyClosedRing{GeomType}(ring::GeomType) where GeomType
        new{GeomType}(ring)
    end
end



function LazyClosedRing(::GI.AbstractCurveTrait, ring::GeomType) where GeomType
    LazyClosedRing{GeomType}(ring)
end

# GeoInterface implementation
GI.geomtrait(::LazyClosedRing) = GI.LinearRingTrait()
GI.is3d(wrapper::LazyClosedRing) = GI.is3d(wrapper.ring)
GI.ismeasured(wrapper::LazyClosedRing) = GI.ismeasured(wrapper.ring)
GI.isclosed(::LazyClosedRing) = true

function GI.npoint(wrapper::LazyClosedRing)
    ring_npoints = GI.npoint(wrapper.ring)
    if GI.getpoint(wrapper.ring, 1) == GI.getpoint(wrapper.ring, ring_npoints)
        return ring_npoints
    else
        return ring_npoints + 1 # account for closing
    end
end

function GI.getpoint(wrapper::LazyClosedRing)
    return LazyClosedRingTuplePointIterator(wrapper)
end

function GI.getpoint(wrapper::LazyClosedRing, i::Integer)
    ring_npoint = GI.npoint(wrapper.ring)
    if i ≤ ring_npoint
        return GI.getpoint(wrapper.ring, i)
    elseif i == ring_npoint + 1
        if GI.getpoint(wrapper.ring, 1) == GI.getpoint(wrapper.ring, ring_npoint)
            return throw(BoundsError(wrapper.ring, i))
        else
            return GI.getpoint(wrapper.ring, 1)
        end
    else
        return throw(BoundsError(wrapper.ring, i))
    end
end

function tuples(wrapper::LazyClosedRing)
    return collect(LazyClosedRingTuplePointIterator(wrapper))
end

struct LazyClosedRingTuplePointIterator{hasZ, hasM, GeomType}
    ring::GeomType
    closed::Bool
end

function LazyClosedRingTuplePointIterator(ring::LazyClosedRing)
    geom = ring.ring
    isclosed = GI.getpoint(geom, 1) == GI.getpoint(geom, GI.npoint(geom))
    return LazyClosedRingTuplePointIterator{GI.is3d(geom), GI.ismeasured(geom), typeof(geom)}(geom, isclosed)
end

# Base iterator interface
Base.IteratorSize(::LazyClosedRingTuplePointIterator) = Base.HasLength()
Base.length(iter::LazyClosedRingTuplePointIterator) = GI.npoint(iter.ring) + iter.closed
Base.IteratorEltype(::LazyClosedRingTuplePointIterator) = Base.HasEltype()
function Base.eltype(::LazyClosedRingTuplePointIterator{hasZ, hasM}) where {hasZ, hasM}
    if !hasZ && !hasM
        Tuple{Float64, Float64}
    elseif hasZ ⊻ hasM
        Tuple{Float64, Float64, Float64}
    else # hasZ && hasM
        Tuple{Float64, Float64, Float64, Float64}
    end
end

function Base.iterate(iter::LazyClosedRingTuplePointIterator)
    return (GI.getpoint(iter.ring, 1), 1)
end

function Base.iterate(iter::LazyClosedRingTuplePointIterator, state)
    ring_npoint = GI.npoint(iter.ring)
    if iter.closed 
        if state == ring_npoint
            return nothing
        else
            return (GI.getpoint(iter.ring, state + 1), state + 1)
        end
    else
        if state < ring_npoint
            return (GI.getpoint(iter.ring, state + 1), state + 1)
        elseif state == ring_npoint 
            return (GI.getpoint(iter.ring, 1), state + 1)
        elseif state == ring_npoint + 1
            return nothing
        else
            throw(BoundsError(iter.ring, state))
        end
    end
end

