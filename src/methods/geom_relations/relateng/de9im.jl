# # DE-9IM matrix, location and dimension codes
export DE9IM
export Mod2Boundary, EndpointBoundary, MultivalentEndpointBoundary, MonovalentEndpointBoundary
#=
Port of code-level concepts from JTS: `Location`, `Dimension`,
`IntersectionMatrix` (pattern matching only), and
`operation/relateng/DimensionLocation.java`.
The `DE9IM` struct is immutable and isbits (`NTuple{9, Int8}`).
=#

# Locations (JTS Location.java)
const LOC_INTERIOR = Int8(0)
const LOC_BOUNDARY = Int8(1)
const LOC_EXTERIOR = Int8(2)
const LOC_NONE     = Int8(-1)

# Dimensions (JTS Dimension.java)
const DIM_FALSE    = Int8(-1)   # 'F'
const DIM_TRUE     = Int8(-2)   # 'T'  (patterns only)
const DIM_DONTCARE = Int8(-3)   # '*'  (patterns only)
const DIM_P = Int8(0)
const DIM_L = Int8(1)
const DIM_A = Int8(2)

function dim_char(d::Integer)
    d == DIM_FALSE    && return 'F'
    d == DIM_TRUE     && return 'T'
    d == DIM_DONTCARE && return '*'
    (0 <= d <= 2)     && return Char('0' + d)
    throw(ArgumentError("invalid dimension code $d"))
end

function dim_code(c::AbstractChar)
    c = uppercase(c)
    c == 'F' && return DIM_FALSE
    c == 'T' && return DIM_TRUE
    c == '*' && return DIM_DONTCARE
    ('0' <= c <= '2') && return Int8(c - '0')
    throw(ArgumentError("invalid DE-9IM character '$c'"))
end

# DimensionLocation codes (JTS DimensionLocation.java, verbatim values)
const DL_EXTERIOR       = LOC_EXTERIOR  # 2
const DL_POINT_INTERIOR = Int8(103)
const DL_LINE_INTERIOR  = Int8(110)
const DL_LINE_BOUNDARY  = Int8(111)
const DL_AREA_INTERIOR  = Int8(120)
const DL_AREA_BOUNDARY  = Int8(121)

dimloc_point(loc::Integer) = loc == LOC_INTERIOR ? DL_POINT_INTERIOR : DL_EXTERIOR
function dimloc_line(loc::Integer)
    loc == LOC_INTERIOR && return DL_LINE_INTERIOR
    loc == LOC_BOUNDARY && return DL_LINE_BOUNDARY
    return DL_EXTERIOR
end
function dimloc_area(loc::Integer)
    loc == LOC_INTERIOR && return DL_AREA_INTERIOR
    loc == LOC_BOUNDARY && return DL_AREA_BOUNDARY
    return DL_EXTERIOR
end
# Explicit switch ports of JTS `DimensionLocation.location` / `.dimension`
# (a `dimloc % 10` shortcut would mis-map `DL_POINT_INTERIOR == 103`).
function dimloc_location(dimloc::Integer)
    (dimloc == DL_POINT_INTERIOR || dimloc == DL_LINE_INTERIOR || dimloc == DL_AREA_INTERIOR) && return LOC_INTERIOR
    (dimloc == DL_LINE_BOUNDARY || dimloc == DL_AREA_BOUNDARY) && return LOC_BOUNDARY
    return LOC_EXTERIOR
end
function dimloc_dimension(dimloc::Integer)
    dimloc == DL_POINT_INTERIOR && return DIM_P
    (dimloc == DL_LINE_INTERIOR || dimloc == DL_LINE_BOUNDARY) && return DIM_L
    (dimloc == DL_AREA_INTERIOR || dimloc == DL_AREA_BOUNDARY) && return DIM_A
    return DIM_FALSE
end
function dimloc_dimension(dimloc::Integer, exterior_dim::Integer)
    dimloc == DL_EXTERIOR && return Int8(exterior_dim)
    return dimloc_dimension(dimloc)
end

"""
    DE9IM

An immutable DE-9IM intersection matrix. Entries are dimension codes
(`DIM_FALSE`, `DIM_P`, `DIM_L`, `DIM_A`) stored row-major over
(Interior, Boundary, Exterior) of A × B, matching the standard string
form `"212101212"`. Construct from a 9-character string or empty
(all-`F`) via `DE9IM()`. Index with `im[locA, locB]`, where the indices are
the JTS location *codes* (`0` = Interior, `1` = Boundary, `2` = Exterior — the
internal `LOC_*` constants), **not** 1-based array positions: `im[0, 0]` is the
Interior/Interior entry.
"""
struct DE9IM
    entries::NTuple{9, Int8}
end
DE9IM() = DE9IM(ntuple(_ -> DIM_FALSE, 9))
function DE9IM(s::AbstractString)
    length(s) == 9 || throw(ArgumentError("DE-9IM string must have 9 characters, got $(repr(s))"))
    return DE9IM(ntuple(i -> dim_code(s[i]), 9))
end

@inline im_index(locA::Integer, locB::Integer) = 3 * Int(locA) + Int(locB) + 1
Base.getindex(im::DE9IM, locA::Integer, locB::Integer) = im.entries[im_index(locA, locB)]
with_entry(im::DE9IM, locA::Integer, locB::Integer, dim::Integer) =
    DE9IM(Base.setindex(im.entries, Int8(dim), im_index(locA, locB)))

# `string(im)` and `"$im"` yield the standard 9-character matrix form via
# this `print`; `show` keeps the constructor form for the REPL.
Base.print(io::IO, im::DE9IM) = join(io, (dim_char(d) for d in im.entries))
Base.show(io::IO, im::DE9IM) = print(io, "DE9IM(\"", string(im), "\")")

"Match a single matrix entry against a pattern code (JTS `IntersectionMatrix.matches`)."
function matches_entry(dim::Int8, pat::Int8)
    pat == DIM_DONTCARE && return true
    pat == DIM_TRUE     && return dim >= 0 || dim == DIM_TRUE
    return dim == pat
end

function im_matches(im::DE9IM, pattern::AbstractString)
    length(pattern) == 9 || throw(ArgumentError("DE-9IM pattern must have 9 characters, got $(repr(pattern))"))
    return all(matches_entry(im.entries[i], dim_code(pattern[i])) for i in 1:9)
end

#=
## Relate queries on a DE-9IM matrix

Ports of the JTS `IntersectionMatrix` named-relationship test methods
(`isContains`, `isWithin`, `isCovers`, `isCoveredBy`, `isCrosses`,
`isEquals`, `isOverlaps`, `isTouches`); `im_matches` above is the port of
`IntersectionMatrix.matches`. These are used as the `value_im`
implementations of the named `IMPredicate` kinds in
`relate_predicates.jl`.
=#

# JTS `IntersectionMatrix.isTrue`: an actual dimension value >= 0 is "true"
# ('T' only occurs in patterns, but is accepted for parity with JTS).
im_is_true(dim::Integer) = dim >= 0 || dim == DIM_TRUE

# Matches `[T*****FF*]`.
is_contains(im::DE9IM) =
    im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
    im[LOC_EXTERIOR, LOC_INTERIOR] == DIM_FALSE &&
    im[LOC_EXTERIOR, LOC_BOUNDARY] == DIM_FALSE

# Matches `[T*F**F***]`.
is_within(im::DE9IM) =
    im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
    im[LOC_INTERIOR, LOC_EXTERIOR] == DIM_FALSE &&
    im[LOC_BOUNDARY, LOC_EXTERIOR] == DIM_FALSE

# Matches `[T*****FF*]`, `[*T****FF*]`, `[***T**FF*]` or `[****T*FF*]`.
function is_covers(im::DE9IM)
    has_point_in_common =
        im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) ||
        im_is_true(im[LOC_INTERIOR, LOC_BOUNDARY]) ||
        im_is_true(im[LOC_BOUNDARY, LOC_INTERIOR]) ||
        im_is_true(im[LOC_BOUNDARY, LOC_BOUNDARY])
    return has_point_in_common &&
        im[LOC_EXTERIOR, LOC_INTERIOR] == DIM_FALSE &&
        im[LOC_EXTERIOR, LOC_BOUNDARY] == DIM_FALSE
end

# Matches `[T*F**F***]`, `[*TF**F***]`, `[**FT*F***]` or `[**F*TF***]`.
function is_coveredby(im::DE9IM)
    has_point_in_common =
        im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) ||
        im_is_true(im[LOC_INTERIOR, LOC_BOUNDARY]) ||
        im_is_true(im[LOC_BOUNDARY, LOC_INTERIOR]) ||
        im_is_true(im[LOC_BOUNDARY, LOC_BOUNDARY])
    return has_point_in_common &&
        im[LOC_INTERIOR, LOC_EXTERIOR] == DIM_FALSE &&
        im[LOC_BOUNDARY, LOC_EXTERIOR] == DIM_FALSE
end

# Matches `[T*T******]` (P/L, P/A, L/A), `[T*****T**]` (L/P, A/P, A/L)
# or `[0********]` (L/L); false for any other dimension combination.
function is_crosses(im::DE9IM, dimA::Integer, dimB::Integer)
    if (dimA == DIM_P && dimB == DIM_L) ||
       (dimA == DIM_P && dimB == DIM_A) ||
       (dimA == DIM_L && dimB == DIM_A)
        return im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
            im_is_true(im[LOC_INTERIOR, LOC_EXTERIOR])
    end
    if (dimA == DIM_L && dimB == DIM_P) ||
       (dimA == DIM_A && dimB == DIM_P) ||
       (dimA == DIM_A && dimB == DIM_L)
        return im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
            im_is_true(im[LOC_EXTERIOR, LOC_INTERIOR])
    end
    if dimA == DIM_L && dimB == DIM_L
        return im[LOC_INTERIOR, LOC_INTERIOR] == DIM_P
    end
    return false
end

# Dimensions equal and matrix matches `[T*F**FFF*]`. (JTS deliberately
# deviates from the SFS pattern `[TFFFTFFFT]` so identical points are equal.)
function is_equals(im::DE9IM, dimA::Integer, dimB::Integer)
    dimA != dimB && return false
    return im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
        im[LOC_INTERIOR, LOC_EXTERIOR] == DIM_FALSE &&
        im[LOC_BOUNDARY, LOC_EXTERIOR] == DIM_FALSE &&
        im[LOC_EXTERIOR, LOC_INTERIOR] == DIM_FALSE &&
        im[LOC_EXTERIOR, LOC_BOUNDARY] == DIM_FALSE
end

# Matches `[T*T***T**]` (P/P, A/A) or `[1*T***T**]` (L/L).
function is_overlaps(im::DE9IM, dimA::Integer, dimB::Integer)
    if (dimA == DIM_P && dimB == DIM_P) || (dimA == DIM_A && dimB == DIM_A)
        return im_is_true(im[LOC_INTERIOR, LOC_INTERIOR]) &&
            im_is_true(im[LOC_INTERIOR, LOC_EXTERIOR]) &&
            im_is_true(im[LOC_EXTERIOR, LOC_INTERIOR])
    end
    if dimA == DIM_L && dimB == DIM_L
        return im[LOC_INTERIOR, LOC_INTERIOR] == DIM_L &&
            im_is_true(im[LOC_INTERIOR, LOC_EXTERIOR]) &&
            im_is_true(im[LOC_EXTERIOR, LOC_INTERIOR])
    end
    return false
end

# Matches `[FT*******]`, `[F**T*****]` or `[F***T****]`;
# false if both geometries are points.
function is_touches(im::DE9IM, dimA::Integer, dimB::Integer)
    if dimA > dimB
        # no need to get transpose because the pattern matrix is symmetrical
        return is_touches(im, dimB, dimA)
    end
    if (dimA == DIM_A && dimB == DIM_A) ||
       (dimA == DIM_L && dimB == DIM_L) ||
       (dimA == DIM_L && dimB == DIM_A) ||
       (dimA == DIM_P && dimB == DIM_A) ||
       (dimA == DIM_P && dimB == DIM_L)
        return im[LOC_INTERIOR, LOC_INTERIOR] == DIM_FALSE &&
            (im_is_true(im[LOC_INTERIOR, LOC_BOUNDARY]) ||
             im_is_true(im[LOC_BOUNDARY, LOC_INTERIOR]) ||
             im_is_true(im[LOC_BOUNDARY, LOC_BOUNDARY]))
    end
    return false
end

# Boundary node rules (JTS BoundaryNodeRule.java). Zero-size structs.
"""
    BoundaryNodeRule

Abstract supertype for rules deciding which endpoints of a linear geometry
are on its boundary, given the number of line ends meeting at the point
(port of JTS `BoundaryNodeRule`). Concrete rules are [`Mod2Boundary`](@ref)
(the OGC SFS default), [`EndpointBoundary`](@ref),
[`MultivalentEndpointBoundary`](@ref), and
[`MonovalentEndpointBoundary`](@ref); each implements
`is_in_boundary(rule, boundary_count)`.
"""
abstract type BoundaryNodeRule end
"""
    Mod2Boundary()

The OGC SFS standard [`BoundaryNodeRule`](@ref) (and the RelateNG default):
an endpoint is on the boundary iff an odd number of line ends meet it
(the "Mod-2 rule"; JTS `Mod2BoundaryNodeRule`).
"""
struct Mod2Boundary <: BoundaryNodeRule end

"""
    EndpointBoundary()

[`BoundaryNodeRule`](@ref) under which any endpoint of a line is on the
boundary, regardless of how many line ends meet there (JTS
`EndPointBoundaryNodeRule`).
"""
struct EndpointBoundary <: BoundaryNodeRule end

"""
    MultivalentEndpointBoundary()

[`BoundaryNodeRule`](@ref) under which an endpoint is on the boundary iff
**more than one** line end meets it (JTS
`MultiValentEndPointBoundaryNodeRule`).
"""
struct MultivalentEndpointBoundary <: BoundaryNodeRule end

"""
    MonovalentEndpointBoundary()

[`BoundaryNodeRule`](@ref) under which an endpoint is on the boundary iff
**exactly one** line end meets it (JTS
`MonoValentEndPointBoundaryNodeRule`).
"""
struct MonovalentEndpointBoundary <: BoundaryNodeRule end

is_in_boundary(::Mod2Boundary, boundary_count::Integer) = isodd(boundary_count)
is_in_boundary(::EndpointBoundary, boundary_count::Integer) = boundary_count > 0
is_in_boundary(::MultivalentEndpointBoundary, boundary_count::Integer) = boundary_count > 1
is_in_boundary(::MonovalentEndpointBoundary, boundary_count::Integer) = boundary_count == 1
