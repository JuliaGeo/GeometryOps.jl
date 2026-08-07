# Minimal FFI binding to s2geography's C ABI (`S2Geography_jll` >= 0.4), used as
# the SPHERICAL oracle for the overlay differential suite.
#
# s2geography wraps Google's S2, so its edges are geodesics — unlike LibGEOS,
# which is planar and therefore cannot validate anything about `Spherical()`.
# This is a test-only binding: nothing here is loaded by GeometryOps itself.
#
# ## The two routes the library offers, and why both are here
#
# 1. **Handle API** (`s2geography_c.h`): `S2GeogCreate` / `S2GeogFactory*` /
#    `S2GeogOp*`. Simple, but PREDICATES ONLY — six op ids, all returning a
#    bool. Used by [`s2_predicate`](@ref).
#
# 2. **Sedona scalar-UDF kernels** (`s2geography/sedona_udf/sedona_extension.h`):
#    the only route that exposes CONSTRUCTIVE overlay (`st_union`,
#    `st_intersection`, `st_difference`, `st_symdifference`) and `st_area`.
#    `S2GeogInitKernels` fills an array of 36 `SedonaCScalarKernel` structs
#    (4 pointers each); each yields a `SedonaCScalarKernelImpl` (5 pointers)
#    that is `init`ed with Arrow argument types and then `execute`d on Arrow
#    arrays. Used by [`s2_overlay`](@ref) and [`s2_area`](@ref).
#
#    nanoarrow is NOT a dependency of this file: `sedona_extension.h` declares
#    `ArrowSchema` / `ArrowArray` inline as plain structs (the Arrow C Data
#    Interface), so Julia builds them directly. Only nanoarrow's *helpers* are
#    static-inline, and we need none of them.
#
# ## The Arrow types the kernels accept
#
# Geography arguments must be `geoarrow.wkb` with SPHERICAL edges. Concretely, a
# binary array (`format = "z"`, int32 offsets) whose schema metadata carries
#
#     ARROW:extension:name     = geoarrow.wkb
#     ARROW:extension:metadata = {"edges":"spherical"}
#
# The `edges` key is load-bearing and is the GeoArrow 0.2 spelling; `edge_type`
# is silently rejected — a kernel that does not recognise its argument types
# reports "does not apply" (a released output schema) rather than an error, so a
# wrong spelling presents as a missing overload, not as a diagnosable failure.
# Results come back as `geoarrow.wkb` (`"z"`) for the four overlay ops and as a
# plain double (`"g"`) for `st_area`.
#
# ## Memory discipline
#
# Every pointer that crosses the boundary sits inside a `GC.@preserve` naming the
# Julia object that owns the bytes, and every `*Create` / `new_impl` / output
# struct has its matching destroy or `release` on both the success and the
# failure path. A mistake here segfaults the process with no stack.

module S2Geog

export s2_available, s2_overlay, s2_area, s2_predicate, S2_RADIUS

const AVAILABLE = if Sys.iswindows()
    # The Windows libs2geography_c loads and answers metadata queries
    # (`s2_versions`, `s2_kernel_names`), but the first Sedona kernel `execute`
    # never returns, hanging CI until the 6-hour job kill. Off until the JLL's
    # Windows build is fixed; oracle coverage remains on the Linux and macOS legs.
    @info "s2geography oracle disabled on Windows (first kernel execution hangs)"
    false
else
    try
        @eval import S2Geography_jll
        true
    catch err
        @info "s2geography oracle unavailable (S2Geography_jll did not load on this platform)" err
        false
    end
end

"""
    s2_available() -> Bool

Whether the s2geography oracle is usable: `S2Geography_jll` loaded and the
platform is not Windows, where kernel execution hangs (see `AVAILABLE`). Suites
gate on this and skip rather than fail.
"""
s2_available() = AVAILABLE

"""
    S2_RADIUS

The sphere radius s2geography measures areas on (`S2Earth::RadiusMeters()`).
GeometryOps uses `GeometryOpsCore.WGS84_EARTH_MEAN_RADIUS = 6371008.8`; the two
differ by 1.2 m, i.e. 3.8e-7 in area. Neither is wrong — divide the ratio out
explicitly rather than widening a tolerance to swallow it.
"""
const S2_RADIUS = 6371010.0

const LIB = AVAILABLE ? S2Geography_jll.libs2geography_c : ""

# ## Arrow C Data Interface (layouts copied from `sedona_extension.h`)

struct CSchema
    format::Ptr{UInt8}; name::Ptr{UInt8}; metadata::Ptr{UInt8}
    flags::Int64; n_children::Int64
    children::Ptr{Ptr{Cvoid}}; dictionary::Ptr{Cvoid}
    release::Ptr{Cvoid}; private_data::Ptr{Cvoid}
end
CSchema() = CSchema(C_NULL, C_NULL, C_NULL, 0, 0, C_NULL, C_NULL, C_NULL, C_NULL)

struct CArray
    length::Int64; null_count::Int64; offset::Int64
    n_buffers::Int64; n_children::Int64
    buffers::Ptr{Ptr{Cvoid}}; children::Ptr{Ptr{Cvoid}}; dictionary::Ptr{Cvoid}
    release::Ptr{Cvoid}; private_data::Ptr{Cvoid}
end
CArray() = CArray(0, 0, 0, 0, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)

# `SedonaCScalarKernel` / `SedonaCScalarKernelImpl`
struct CKernel
    function_name::Ptr{Cvoid}; new_impl::Ptr{Cvoid}
    release::Ptr{Cvoid}; private_data::Ptr{Cvoid}
end
struct CImpl
    init::Ptr{Cvoid}; execute::Ptr{Cvoid}; get_last_error::Ptr{Cvoid}
    release::Ptr{Cvoid}; private_data::Ptr{Cvoid}
end
CImpl() = CImpl(C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)

# ## Arrow input plumbing
#
# Input arrays and schemas are owned by Julia and outlive every call that sees
# them, so their release callback is a no-op rather than NULL (a NULL release
# marks a struct as already released, which the consumer is entitled to reject).

_noop_release(::Ptr{Cvoid}) = nothing
const NOOP_RELEASE = @cfunction(_noop_release, Cvoid, (Ptr{Cvoid},))

function _encode_metadata(pairs)
    io = IOBuffer()
    write(io, Int32(length(pairs)))
    for (k, v) in pairs
        write(io, Int32(sizeof(k))); write(io, codeunits(k))
        write(io, Int32(sizeof(v))); write(io, codeunits(v))
    end
    return take!(io)
end

const FMT_BINARY = Vector{UInt8}("z\0")
const EMPTY_NAME = Vector{UInt8}("\0")
const WKB_METADATA = _encode_metadata([
    "ARROW:extension:name" => "geoarrow.wkb",
    "ARROW:extension:metadata" => "{\"edges\":\"spherical\"}",
])

_wkb_schema() = CSchema(pointer(FMT_BINARY), pointer(EMPTY_NAME), pointer(WKB_METADATA),
                        2 #= ARROW_FLAG_NULLABLE =#, 0, C_NULL, C_NULL, NOOP_RELEASE, C_NULL)

# A length-1 binary array viewing `wkb` with no copy. `offs` and `bufs` are
# caller-owned scratch that must outlive the call.
function _wkb_array(wkb::Vector{UInt8}, offs::Vector{Int32}, bufs::Vector{Ptr{Cvoid}})
    offs[1] = Int32(0); offs[2] = Int32(length(wkb))
    bufs[1] = C_NULL                            # validity: no nulls
    bufs[2] = Ptr{Cvoid}(pointer(offs))
    bufs[3] = Ptr{Cvoid}(pointer(wkb))
    return CArray(1, 0, 0, 3, 0, pointer(bufs), C_NULL, C_NULL, NOOP_RELEASE, C_NULL)
end

_release!(r::Base.RefValue{CSchema}) = r[].release == C_NULL ? nothing :
    ccall(r[].release, Cvoid, (Ptr{CSchema},), Base.unsafe_convert(Ptr{CSchema}, r))
_release!(r::Base.RefValue{CArray}) = r[].release == C_NULL ? nothing :
    ccall(r[].release, Cvoid, (Ptr{CArray},), Base.unsafe_convert(Ptr{CArray}, r))
_release!(r::Base.RefValue{CImpl}) = r[].release == C_NULL ? nothing :
    ccall(r[].release, Cvoid, (Ptr{CImpl},), Base.unsafe_convert(Ptr{CImpl}, r))

# Copy element 1 out of a result array. The bytes belong to the kernel, so they
# are copied out before the array is released.
function _read_element(arr::CArray, fmt::String)
    arr.length == 1 || error("expected a 1-row result, got $(arr.length)")
    arr.null_count == 0 || error("s2geography returned a null result")
    bufs = unsafe_wrap(Array, Ptr{Ptr{Cvoid}}(arr.buffers), Int(arr.n_buffers))
    i = 1 + arr.offset
    if fmt == "g"
        return unsafe_load(Ptr{Float64}(bufs[2]), i)
    elseif fmt == "z" || fmt == "Z"
        T = fmt == "z" ? Int32 : Int64
        lo = unsafe_load(Ptr{T}(bufs[2]), i)
        hi = unsafe_load(Ptr{T}(bufs[2]), i + 1)
        n = Int(hi - lo)
        out = Vector{UInt8}(undef, n)
        GC.@preserve out unsafe_copyto!(pointer(out), Ptr{UInt8}(bufs[3]) + lo, n)
        return out
    end
    error("unhandled Arrow result format \"$fmt\"")
end

# ## Kernel registry (filled by `_init_kernels!` at the bottom, once)

const NKERNEL = Ref(0)
const KERNELS = Ref(CKernel[])
const KERNEL_INDEX = Dict{String, Vector{Int}}()

# `function_name` reads `self->private_data`, so a kernel must be addressed where
# it lives in the array, never through a copy.
function _kernel_name(i)
    ks = KERNELS[]
    GC.@preserve ks unsafe_string(ccall(ks[i].function_name, Cstring, (Ptr{Cvoid},), pointer(ks, i)))
end

function _init_kernels!()
    n = Int(ccall((:S2GeogNumKernels, LIB), Csize_t, ()))
    ks = Vector{CKernel}(undef, n)
    rc = GC.@preserve ks ccall((:S2GeogInitKernels, LIB), Cint,
        (Ptr{Cvoid}, Csize_t, Cint), pointer(ks), sizeof(ks),
        1)                                   # S2GEOGRAPHY_KERNEL_FORMAT_SEDONA_UDF
    rc == 0 || error("S2GeogInitKernels failed with $rc")
    NKERNEL[] = n
    KERNELS[] = ks
    empty!(KERNEL_INDEX)
    for i in 1:n
        push!(get!(Vector{Int}, KERNEL_INDEX, _kernel_name(i)), i)
    end
    return nothing
end

_last_error(impl, pimpl) =
    unsafe_string(ccall(impl[].get_last_error, Cstring, (Ptr{CImpl},), pimpl))

"""
    _run_kernel(name, wkbs) -> Vector{UInt8} | Float64

Run the scalar UDF `name` over one row of `length(wkbs)` geography arguments.
Each registered overload of `name` is tried in turn; an overload that reports
"does not apply" (a released output schema) is skipped.
"""
function _run_kernel(name::String, wkbs::Vector{Vector{UInt8}})
    idxs = get(KERNEL_INDEX, name, Int[])
    isempty(idxs) && error("s2geography exports no kernel named $name")
    kernels = KERNELS[]
    n = length(wkbs)
    schemas = [Ref(_wkb_schema()) for _ in 1:n]
    offsets = [Int32[0, 0] for _ in 1:n]
    bufptrs = [Vector{Ptr{Cvoid}}(undef, 3) for _ in 1:n]
    GC.@preserve wkbs schemas offsets bufptrs kernels begin
        arrays = [Ref(_wkb_array(wkbs[i], offsets[i], bufptrs[i])) for i in 1:n]
        argtypes = Ptr{Cvoid}[Base.unsafe_convert(Ptr{CSchema}, s) for s in schemas]
        args = Ptr{Cvoid}[Base.unsafe_convert(Ptr{CArray}, a) for a in arrays]
        GC.@preserve arrays argtypes args begin
            for i in idxs
                impl = Ref(CImpl()); outs = Ref(CSchema()); outa = Ref(CArray())
                GC.@preserve impl outs outa begin
                    pimpl = Base.unsafe_convert(Ptr{CImpl}, impl)
                    pouts = Base.unsafe_convert(Ptr{CSchema}, outs)
                    pouta = Base.unsafe_convert(Ptr{CArray}, outa)
                    ccall(kernels[i].new_impl, Cvoid, (Ptr{Cvoid}, Ptr{CImpl}),
                          pointer(kernels, i), pimpl)
                    try
                        rc = ccall(impl[].init, Cint,
                            (Ptr{CImpl}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Int64, Ptr{CSchema}),
                            pimpl, pointer(argtypes), C_NULL, Int64(n), pouts)
                        rc == 0 || error("$name init failed: " * _last_error(impl, pimpl))
                        #-- a released output schema means "this overload does not
                        #-- apply to these argument types"
                        outs[].release == C_NULL && continue
                        fmt = unsafe_string(outs[].format)
                        rc = ccall(impl[].execute, Cint,
                            (Ptr{CImpl}, Ptr{Ptr{Cvoid}}, Int64, Int64, Ptr{CArray}),
                            pimpl, pointer(args), Int64(n), Int64(1), pouta)
                        rc == 0 || error("$name failed: " * _last_error(impl, pimpl))
                        return _read_element(outa[], fmt)
                    finally
                        _release!(outa); _release!(outs); _release!(impl)
                    end
                end
            end
        end
    end
    error("no overload of $name accepted geoarrow.wkb/spherical arguments")
end

# ## Public surface — constructive ops (route 2)

const OVERLAY_KERNEL = Dict(
    :intersection => "st_intersection", :union => "st_union",
    :difference => "st_difference", :symdifference => "st_symdifference")

"""
    s2_overlay(op, wkb_a, wkb_b) -> Vector{UInt8}

s2geography's `op` (`:intersection`, `:union`, `:difference`, `:symdifference`)
over two WKB buffers, as a WKB buffer. Edges are geodesics on both sides.
"""
function s2_overlay(op::Symbol, wkb_a::Vector{UInt8}, wkb_b::Vector{UInt8})
    k = get(OVERLAY_KERNEL, op) do
        error("unknown overlay op $op")
    end
    return _run_kernel(k, Vector{UInt8}[wkb_a, wkb_b])::Vector{UInt8}
end

"""
    s2_area(wkb) -> Float64

Geodesic area in m² on a sphere of radius [`S2_RADIUS`](@ref).
"""
s2_area(wkb::Vector{UInt8}) = _run_kernel("st_area", Vector{UInt8}[wkb])::Float64

# ## Public surface — predicates (the handle API, route 1)

const PREDICATE_OP = Dict(:intersects => 1, :contains => 2, :within => 3,
                          :equals => 4, :disjoint => 6)

_error_message(err) =
    unsafe_string(ccall((:S2GeogErrorGetMessage, LIB), Cstring, (Ptr{Cvoid},), err))

"""
    s2_predicate(op, wkb_a, wkb_b) -> Bool

One of s2geography's handle-API predicates (`:intersects`, `:contains`,
`:within`, `:equals`, `:disjoint`) over two WKB buffers. `DISTANCE_WITHIN` is
omitted: it needs the three-argument evaluator and a threshold.
"""
function s2_predicate(op::Symbol, wkb_a::Vector{UInt8}, wkb_b::Vector{UInt8})
    id = get(PREDICATE_OP, op) do
        error("unknown predicate $op")
    end
    err = Ref(Ptr{Cvoid}(C_NULL)); factory = Ref(Ptr{Cvoid}(C_NULL))
    ga = Ref(Ptr{Cvoid}(C_NULL)); gb = Ref(Ptr{Cvoid}(C_NULL)); opp = Ref(Ptr{Cvoid}(C_NULL))
    ccall((:S2GeogErrorCreate, LIB), Cint, (Ptr{Ptr{Cvoid}},), err) == 0 ||
        error("S2GeogErrorCreate failed")
    try
        ccall((:S2GeogFactoryCreate, LIB), Cint, (Ptr{Ptr{Cvoid}},), factory) == 0 ||
            error("S2GeogFactoryCreate failed")
        ccall((:S2GeogCreate, LIB), Cint, (Ptr{Ptr{Cvoid}},), ga) == 0 || error("S2GeogCreate failed")
        ccall((:S2GeogCreate, LIB), Cint, (Ptr{Ptr{Cvoid}},), gb) == 0 || error("S2GeogCreate failed")
        ccall((:S2GeogOpCreate, LIB), Cint, (Ptr{Ptr{Cvoid}}, Cint), opp, id) == 0 ||
            error("S2GeogOpCreate failed")
        #-- InitFromWkbNonOwning does not copy, so both buffers must stay pinned
        #-- until the op that reads the geographies has returned
        GC.@preserve wkb_a wkb_b begin
            for (buf, g) in ((wkb_a, ga), (wkb_b, gb))
                rc = ccall((:S2GeogFactoryInitFromWkbNonOwning, LIB), Cint,
                    (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{Cvoid}, Ptr{Cvoid}),
                    factory[], pointer(buf), length(buf), g[], err[])
                rc == 0 || error("WKB import failed: " * _error_message(err[]))
            end
            rc = ccall((:S2GeogOpEvalGeogGeog, LIB), Cint,
                (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), opp[], ga[], gb[], err[])
            rc == 0 || error("$op failed: " * _error_message(err[]))
        end
        return ccall((:S2GeogOpGetInt, LIB), Int64, (Ptr{Cvoid},), opp[]) != 0
    finally
        opp[] == C_NULL || ccall((:S2GeogOpDestroy, LIB), Cvoid, (Ptr{Cvoid},), opp[])
        ga[] == C_NULL || ccall((:S2GeogDestroy, LIB), Cvoid, (Ptr{Cvoid},), ga[])
        gb[] == C_NULL || ccall((:S2GeogDestroy, LIB), Cvoid, (Ptr{Cvoid},), gb[])
        factory[] == C_NULL || ccall((:S2GeogFactoryDestroy, LIB), Cvoid, (Ptr{Cvoid},), factory[])
        ccall((:S2GeogErrorDestroy, LIB), Cvoid, (Ptr{Cvoid},), err[])
    end
end

"""
    s2_versions() -> NamedTuple

Versions of the libraries s2geography was built against, for the run banner.
"""
s2_versions() = (
    s2geometry = unsafe_string(ccall((:S2GeogS2GeometryVersion, LIB), Cstring, ())),
    geoarrow = unsafe_string(ccall((:S2GeogGeoArrowVersion, LIB), Cstring, ())),
    nanoarrow = unsafe_string(ccall((:S2GeogNanoarrowVersion, LIB), Cstring, ())),
)

"""
    s2_kernel_names() -> Vector{String}

Names of the exported scalar UDFs, in registry order.
"""
s2_kernel_names() = [_kernel_name(i) for i in 1:NKERNEL[]]

# A top-level statement, so every method above is already in the caller's world.
AVAILABLE && _init_kernels!()

end  # module
