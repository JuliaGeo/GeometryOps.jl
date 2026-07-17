# # Fast WKB codec for polygonal geometry
export parse_wkb, write_wkb

#=
## What is this?

A dependency-light, single-pass reader and writer for **Well-Known Binary (WKB)**
restricted to polygonal geometry (`Polygon` and `MultiPolygon`).  It exists to move
geometry in and out of C libraries (s2geography's C API, LibGEOS) as `geoarrow.wkb`
(standard ISO WKB) during differential-validation testing, without pulling in a full
WKB stack.

* [`parse_wkb`](@ref) reads a byte buffer and returns a `GI.Polygon` or `GI.MultiPolygon`
  whose rings are `Vector{Tuple{Float64,Float64}}` — the same nested-tuple representation
  `GO.tuples` produces, so results compare `==` coordinate-exact against it.
* [`write_wkb`](@ref) serializes any GeoInterface `PolygonTrait`/`MultiPolygonTrait`
  geometry to little-endian ISO WKB.

Parsing is the hot path (results come back from C as WKB), so the reader is a single
pass over the buffer using pointer loads — no `IOBuffer`, no per-value `read`, no
intermediate geometry objects.  Every count is validated against the remaining buffer
length *before* any allocation, so a truncated or hostile buffer throws a clean error
instead of segfaulting or OOMing on a giant preallocation.

## WKB layout (ISO)

Each geometry is: a **byte-order** byte (`0x00` big-endian/XDR, `0x01` little-endian/NDR),
a `UInt32` **type code**, then a type-specific payload.  A `MultiPolygon` payload is a
`UInt32` count followed by full sub-geometries — **each carrying its own byte-order byte
and type code**, so mixed endianness within one buffer is legal and handled.

Dimensionality is detected two ways and rejected (2D only): ISO codes (`type + 1000` Z,
`+ 2000` M, `+ 3000` ZM) and EWKB high-bit flags (`0x80000000` Z, `0x40000000` M).  The
EWKB SRID flag (`0x20000000`) is accepted and its 4-byte SRID skipped.
=#

import GeoInterface as GI

# Concrete geometry types.  These match exactly what `GO.tuples` returns for 2D
# polygonal geometry, so `parse_wkb(...) == GO.tuples(geom)` holds structurally.
# Building them through the inner constructors (rather than the public `GI.Polygon`
# etc.) skips validation and, crucially, works for empty rings/polygons/multipolygons
# (the public constructors call `first(geom)` and throw on empty input).
const _WKBRing  = Vector{Tuple{Float64,Float64}}
const _WKBLinRing = GI.LinearRing{false,false,_WKBRing,Nothing,Nothing}
const _WKBRings = Vector{_WKBLinRing}
const _WKBPoly  = GI.Polygon{false,false,_WKBRings,Nothing,Nothing}
const _WKBPolys = Vector{_WKBPoly}
const _WKBMPoly = GI.MultiPolygon{false,false,_WKBPolys,Nothing,Nothing}

@inline _wkb_linearring(coords::_WKBRing) = _WKBLinRing(coords, nothing, nothing)
@inline _wkb_polygon(rings::_WKBRings) = _WKBPoly(rings, nothing, nothing)
@inline _wkb_multipolygon(polys::_WKBPolys) = _WKBMPoly(polys, nothing, nothing)

# WKB integer type codes -> human names, for error messages.
function _wkb_type_name(code::Integer)
    code == 1 ? "Point" :
    code == 2 ? "LineString" :
    code == 3 ? "Polygon" :
    code == 4 ? "MultiPoint" :
    code == 5 ? "MultiLineString" :
    code == 6 ? "MultiPolygon" :
    code == 7 ? "GeometryCollection" :
    "unknown"
end

# Cold error paths, kept out of the hot loop.
@noinline _wkb_truncated(need, have) =
    throw(ArgumentError("Truncated WKB buffer: need $(need) more byte(s) but only $(have) remain"))
@noinline _wkb_bad_order(b) =
    throw(ArgumentError("Invalid WKB byte-order flag $(Int(b)); expected 0 (big-endian) or 1 (little-endian)"))
@noinline _wkb_zm(rawtype) =
    throw(ArgumentError("parse_wkb supports 2D geometries only, but WKB type code $(Int(rawtype)) declares Z and/or M coordinates"))
@noinline _wkb_unsupported(code) =
    throw(ArgumentError("parse_wkb only supports Polygon (3) and MultiPolygon (6); found WKB type $(code) ($(_wkb_type_name(code)))"))
@noinline _wkb_mp_child(code) =
    throw(ArgumentError("A MultiPolygon element must be a Polygon (3), but found WKB type $(code) ($(_wkb_type_name(code)))"))

# ## Low-level reads.
# `p0` points at byte 1; `pos` is a 0-based byte offset; `n` is the buffer length.
# `ltoh`/`ntoh` make endianness host-independent: little-endian data goes through
# `ltoh`, big-endian through `ntoh`, each a no-op or `bswap` depending on host.

@inline function _read_u8(p0::Ptr{UInt8}, pos::Int, n::Int)
    pos + 1 <= n || _wkb_truncated(1, n - pos)
    return unsafe_load(p0 + pos), pos + 1
end

@inline function _read_u32(p0::Ptr{UInt8}, pos::Int, n::Int, le::Bool)
    pos + 4 <= n || _wkb_truncated(4, n - pos)
    v = unsafe_load(Ptr{UInt32}(p0 + pos))
    return (le ? ltoh(v) : ntoh(v)), pos + 4
end

# No bounds check here: callers validate the whole run of coordinates up front.
@inline function _read_f64(p0::Ptr{UInt8}, pos::Int, le::Bool)
    u = unsafe_load(Ptr{UInt64}(p0 + pos))
    return reinterpret(Float64, le ? ltoh(u) : ntoh(u)), pos + 8
end

# Read a geometry header (byte order + type + optional SRID); return `(geomcode, le, pos)`.
@inline function _read_header(p0::Ptr{UInt8}, pos::Int, n::Int)
    order, pos = _read_u8(p0, pos, n)
    le = order == 0x01 ? true : order == 0x00 ? false : _wkb_bad_order(order)
    rawtype, pos = _read_u32(p0, pos, n, le)
    has_z    = (rawtype & 0x80000000) != 0
    has_m    = (rawtype & 0x40000000) != 0
    has_srid = (rawtype & 0x20000000) != 0
    base = rawtype & 0x1fffffff          # strip the three EWKB flag bits
    isodim  = base ÷ UInt32(1000)        # ISO Z/M encoding lives in the thousands digit
    geomcode = base % UInt32(1000)
    (has_z || has_m || isodim != 0) && _wkb_zm(rawtype)
    if has_srid
        _, pos = _read_u32(p0, pos, n, le)   # accept and discard the 4-byte SRID
    end
    return Int(geomcode), le, pos
end

# Parse a Polygon body (header already consumed); return `(polygon, pos)`.
@inline function _read_polygon_body(p0::Ptr{UInt8}, pos::Int, n::Int, le::Bool)
    nrings32, pos = _read_u32(p0, pos, n, le)
    nrings = Int(nrings32)
    # Each ring is at least its own 4-byte point count; bound the allocation first.
    Int64(pos) + Int64(nrings) * 4 <= n || _wkb_truncated(nrings * 4, n - pos)
    rings = Vector{_WKBLinRing}(undef, nrings)
    @inbounds for r in 1:nrings
        npts32, pos = _read_u32(p0, pos, n, le)
        npts = Int(npts32)
        need = Int64(npts) * 16
        Int64(pos) + need <= n || _wkb_truncated(need, n - pos)
        coords = Vector{Tuple{Float64,Float64}}(undef, npts)
        for i in 1:npts
            x, pos = _read_f64(p0, pos, le)
            y, pos = _read_f64(p0, pos, le)
            coords[i] = (x, y)
        end
        rings[r] = _wkb_linearring(coords)
    end
    return _wkb_polygon(rings), pos
end

# Parse a MultiPolygon body (header already consumed); return `(multipolygon, pos)`.
@inline function _read_multipolygon_body(p0::Ptr{UInt8}, pos::Int, n::Int, le::Bool)
    ngeoms32, pos = _read_u32(p0, pos, n, le)
    ngeoms = Int(ngeoms32)
    # Each sub-polygon is at least order(1) + type(4) + nrings(4) = 9 bytes.
    Int64(pos) + Int64(ngeoms) * 9 <= n || _wkb_truncated(ngeoms * 9, n - pos)
    polys = Vector{_WKBPoly}(undef, ngeoms)
    @inbounds for g in 1:ngeoms
        geomcode, le2, pos = _read_header(p0, pos, n)   # per-element byte order + type
        geomcode == 3 || _wkb_mp_child(geomcode)
        poly, pos = _read_polygon_body(p0, pos, n, le2)
        polys[g] = poly
    end
    return _wkb_multipolygon(polys), pos
end

"""
    parse_wkb(bytes::AbstractVector{UInt8})

Parse ISO/EWKB Well-Known Binary into a `GI.Polygon` or `GI.MultiPolygon` whose rings
are `Vector{Tuple{Float64,Float64}}` (the `GO.tuples` representation).

Only 2D `Polygon` (type 3) and `MultiPolygon` (type 6) are supported, including empty
ones.  Big-endian and little-endian buffers are both handled, as is mixed endianness
across the sub-geometries of a `MultiPolygon`.  The EWKB SRID flag is accepted and the
SRID ignored; any Z/M dimension (ISO or EWKB flavor) throws.  Malformed or truncated
buffers throw an `ArgumentError` rather than crashing.
"""
function parse_wkb(bytes::AbstractVector{UInt8})
    n = length(bytes)
    GC.@preserve bytes begin
        p0 = pointer(bytes)
        geomcode, le, pos = _read_header(p0, 0, n)
        if geomcode == 3
            poly, _ = _read_polygon_body(p0, pos, n, le)
            return poly
        elseif geomcode == 6
            mp, _ = _read_multipolygon_body(p0, pos, n, le)
            return mp
        else
            _wkb_unsupported(geomcode)
        end
    end
end

# ## Writing.
# Little-endian ISO WKB.  `htol` maps host order to little-endian (a no-op or `bswap`).

@inline function _put_u8!(p0::Ptr{UInt8}, pos::Int, v::UInt8)
    unsafe_store!(p0 + pos, v)
    return pos + 1
end
@inline function _put_u32!(p0::Ptr{UInt8}, pos::Int, v::UInt32)
    unsafe_store!(Ptr{UInt32}(p0 + pos), htol(v))
    return pos + 4
end
@inline function _put_f64!(p0::Ptr{UInt8}, pos::Int, v::Float64)
    unsafe_store!(Ptr{UInt64}(p0 + pos), htol(reinterpret(UInt64, v)))
    return pos + 8
end

# Byte size of a Polygon / MultiPolygon in little-endian ISO WKB.
function _polygon_wkb_size(poly)
    sz = 1 + 4 + 4                       # order + type + ring count
    for ring in GI.getgeom(poly)
        sz += 4 + Int(GI.npoint(ring)) * 16
    end
    return sz
end
function _multipolygon_wkb_size(mp)
    sz = 1 + 4 + 4                       # order + type + polygon count
    for poly in GI.getgeom(mp)
        sz += _polygon_wkb_size(poly)
    end
    return sz
end

function _write_polygon!(p0::Ptr{UInt8}, pos::Int, poly)
    pos = _put_u8!(p0, pos, 0x01)                     # little-endian
    pos = _put_u32!(p0, pos, UInt32(3))               # Polygon
    pos = _put_u32!(p0, pos, UInt32(GI.ngeom(poly)))
    for ring in GI.getgeom(poly)
        pos = _put_u32!(p0, pos, UInt32(GI.npoint(ring)))
        for pt in GI.getpoint(ring)
            pos = _put_f64!(p0, pos, Float64(GI.x(pt)))
            pos = _put_f64!(p0, pos, Float64(GI.y(pt)))
        end
    end
    return pos
end

function _write_multipolygon!(p0::Ptr{UInt8}, pos::Int, mp)
    pos = _put_u8!(p0, pos, 0x01)                     # little-endian
    pos = _put_u32!(p0, pos, UInt32(6))               # MultiPolygon
    pos = _put_u32!(p0, pos, UInt32(GI.ngeom(mp)))
    for poly in GI.getgeom(mp)
        pos = _write_polygon!(p0, pos, poly)
    end
    return pos
end

"""
    write_wkb(geom)::Vector{UInt8}

Serialize a 2D GeoInterface `PolygonTrait` or `MultiPolygonTrait` geometry to
little-endian ISO WKB (type 3 / 6).  Doubles are written bit-for-bit; the output buffer
is preallocated from an exact size computation.  Throws `ArgumentError` for any other
geometry.
"""
function write_wkb(geom)
    t = GI.trait(geom)
    if t isa GI.PolygonTrait
        out = Vector{UInt8}(undef, _polygon_wkb_size(geom))
        GC.@preserve out _write_polygon!(pointer(out), 0, geom)
        return out
    elseif t isa GI.MultiPolygonTrait
        out = Vector{UInt8}(undef, _multipolygon_wkb_size(geom))
        GC.@preserve out _write_multipolygon!(pointer(out), 0, geom)
        return out
    else
        throw(ArgumentError("write_wkb only supports Polygon and MultiPolygon geometries, got trait $(t)"))
    end
end
