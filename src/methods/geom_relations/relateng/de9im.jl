# # DE-9IM matrix, location and dimension codes
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
(all-`F`) via `DE9IM()`. Index with `im[locA, locB]`.
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

Base.string(im::DE9IM) = join(dim_char(d) for d in im.entries)
Base.show(io::IO, im::DE9IM) = print(io, "DE9IM(\"", string(im), "\")")

"Match a single matrix entry against a pattern code (JTS `IntersectionMatrix.matches`)."
function matches_entry(dim::Int8, pat::Int8)
    pat == DIM_DONTCARE && return true
    pat == DIM_TRUE     && return dim >= 0
    return dim == pat
end

function matches(im::DE9IM, pattern::AbstractString)
    length(pattern) == 9 || throw(ArgumentError("DE-9IM pattern must have 9 characters, got $(repr(pattern))"))
    return all(matches_entry(im.entries[i], dim_code(pattern[i])) for i in 1:9)
end
