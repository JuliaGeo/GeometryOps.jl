# # Fast WKB codec for polygonal geometry
export parse_wkb, write_wkb, WKBParseError, WKBWriteError

#=
## What is this?

A single-pass codec for 2D polygonal **Well-Known Binary**.  [`parse_wkb`](@ref) returns
a `GI.Polygon`/`GI.MultiPolygon` whose rings are `Vector{Tuple{Float64,Float64}}`, equal
to `GO.tuples` output; [`write_wkb`](@ref) emits little-endian ISO WKB.

## WKB layout (ISO)

A geometry is a **byte-order** byte (`0x00` big-endian, `0x01` little-endian), a `UInt32`
**type code**, then a payload.  A `MultiPolygon` payload is a `UInt32` count followed by
whole sub-geometries, each with its own byte-order byte and type code.  Z/M dimensions are
encoded either in the ISO type code (`type + 1000` Z, `+ 2000` M, `+ 3000` ZM) or in EWKB
flag bits (`0x80000000` Z, `0x40000000` M), and are rejected here; the EWKB SRID flag
(`0x20000000`) is accepted and its 4-byte payload skipped.
=#

import GeoInterface as GI

# Concrete types matching `GO.tuples` output.  Their inner constructors skip validation
# and accept empty rings, polygons and multipolygons.
const _WKBLinearRing = GI.LinearRing{false,false,Vector{Tuple{Float64,Float64}},Nothing,Nothing}
const _WKBPolygon = GI.Polygon{false,false,Vector{_WKBLinearRing},Nothing,Nothing}
const _WKBMultiPolygon = GI.MultiPolygon{false,false,Vector{_WKBPolygon},Nothing,Nothing}

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

"""
    abstract type WKBParseError <: Exception

Supertype of every error [`parse_wkb`](@ref) throws on malformed input: `WKBTruncatedError`,
`WKBByteOrderError`, `WKBDimensionError`, `WKBUnsupportedTypeError`, `WKBMultiPolygonElementError`.
"""
abstract type WKBParseError <: Exception end

"""
    WKBWriteError <: Exception

Error thrown by [`write_wkb`](@ref) for geometry it cannot serialize, i.e. anything that
is not a 2D `Polygon` or `MultiPolygon`.
"""
struct WKBWriteError <: Exception
    msg::String
end
Base.showerror(io::IO, e::WKBWriteError) = print(io, "WKBWriteError: ", e.msg)

# A field runs past the end of the buffer.
struct WKBTruncatedError <: WKBParseError
    needed::Int
    available::Int
end
Base.showerror(io::IO, e::WKBTruncatedError) = print(io,
    "WKBTruncatedError: need $(e.needed) more byte(s) but only $(e.available) remain")

# The byte-order flag is neither 0x00 (big-endian) nor 0x01 (little-endian).
struct WKBByteOrderError <: WKBParseError
    flag::UInt8
end
Base.showerror(io::IO, e::WKBByteOrderError) = print(io,
    "WKBByteOrderError: byte-order flag $(Int(e.flag)) is neither 0 (big-endian) nor 1 (little-endian)")

# The type code declares Z and/or M coordinates, in either the ISO or the EWKB encoding.
struct WKBDimensionError <: WKBParseError
    rawtype::UInt32
end
Base.showerror(io::IO, e::WKBDimensionError) = print(io,
    "WKBDimensionError: type code $(Int(e.rawtype)) declares Z and/or M coordinates, but only 2D is supported")

# The top-level geometry is neither a Polygon nor a MultiPolygon.
struct WKBUnsupportedTypeError <: WKBParseError
    code::Int
end
Base.showerror(io::IO, e::WKBUnsupportedTypeError) = print(io,
    "WKBUnsupportedTypeError: type $(e.code) ($(_wkb_type_name(e.code))) is not Polygon (3) or MultiPolygon (6)")

# A MultiPolygon holds a sub-geometry that is not a Polygon.
struct WKBMultiPolygonElementError <: WKBParseError
    code::Int
end
Base.showerror(io::IO, e::WKBMultiPolygonElementError) = print(io,
    "WKBMultiPolygonElementError: a MultiPolygon element must be a Polygon (3), but found type " *
    "$(e.code) ($(_wkb_type_name(e.code)))")

# Cold throw paths, kept out of the inlined readers.
@noinline _wkb_truncated(need, have) = throw(WKBTruncatedError(need, have))
@noinline _wkb_bad_order(b) = throw(WKBByteOrderError(b))
@noinline _wkb_zm(rawtype) = throw(WKBDimensionError(rawtype))
@noinline _wkb_unsupported(code) = throw(WKBUnsupportedTypeError(code))
@noinline _wkb_mp_child(code) = throw(WKBMultiPolygonElementError(code))

# ## Low-level reads
# `p0` points at byte 1, `pos` is a 0-based offset from it, `n` is the buffer length.
# `ltoh`/`ntoh` decode little-/big-endian data on either host.

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
    rings = Vector{_WKBLinearRing}(undef, nrings)
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
        rings[r] = _WKBLinearRing(coords, nothing, nothing)
    end
    return _WKBPolygon(rings, nothing, nothing), pos
end

# Parse a MultiPolygon body (header already consumed); return `(multipolygon, pos)`.
@inline function _read_multipolygon_body(p0::Ptr{UInt8}, pos::Int, n::Int, le::Bool)
    ngeoms32, pos = _read_u32(p0, pos, n, le)
    ngeoms = Int(ngeoms32)
    # Each sub-polygon is at least order(1) + type(4) + nrings(4) = 9 bytes.
    Int64(pos) + Int64(ngeoms) * 9 <= n || _wkb_truncated(ngeoms * 9, n - pos)
    polys = Vector{_WKBPolygon}(undef, ngeoms)
    @inbounds for g in 1:ngeoms
        geomcode, le2, pos = _read_header(p0, pos, n)   # per-element byte order + type
        geomcode == 3 || _wkb_mp_child(geomcode)
        poly, pos = _read_polygon_body(p0, pos, n, le2)
        polys[g] = poly
    end
    return _WKBMultiPolygon(polys, nothing, nothing), pos
end

"""
    parse_wkb(bytes::AbstractVector{UInt8})

Parse ISO/EWKB Well-Known Binary into a `GI.Polygon` or `GI.MultiPolygon` whose rings
are `Vector{Tuple{Float64,Float64}}` (the `GO.tuples` representation).

Only 2D `Polygon` (type 3) and `MultiPolygon` (type 6) are supported, including empty
ones.  Big-endian and little-endian buffers are both handled, as is mixed endianness
across the sub-geometries of a `MultiPolygon`.  The EWKB SRID flag is accepted and the
SRID ignored; any Z/M dimension (ISO or EWKB flavor) throws.  Malformed or truncated
buffers throw a [`WKBParseError`](@ref) rather than crashing.
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

# ## Writing
# Little-endian ISO WKB, emitted into an exactly-sized preallocated buffer.

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
is preallocated from an exact size computation.  Any other geometry, and any geometry
carrying Z or M coordinates, throws a [`WKBWriteError`](@ref).
"""
function write_wkb(geom)
    t = GI.trait(geom)
    (t isa GI.PolygonTrait || t isa GI.MultiPolygonTrait) ||
        throw(WKBWriteError("only Polygon and MultiPolygon geometries can be written, got trait $(t)"))
    (GI.is3d(geom) || GI.ismeasured(geom)) &&
        throw(WKBWriteError("only 2D geometries can be written, but this geometry has Z and/or M coordinates"))
    if t isa GI.PolygonTrait
        out = Vector{UInt8}(undef, _polygon_wkb_size(geom))
        GC.@preserve out _write_polygon!(pointer(out), 0, geom)
    else
        out = Vector{UInt8}(undef, _multipolygon_wkb_size(geom))
        GC.@preserve out _write_multipolygon!(pointer(out), 0, geom)
    end
    return out
end
