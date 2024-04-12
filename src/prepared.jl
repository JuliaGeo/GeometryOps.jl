struct Prepared{Pa,Pr}
    parent::Pa
    preparations::Pr
end

Base.parent(x::Prepared) = x.parent
@inline getprep(p::Prepared, x::Symbol) = getpropery(p.preparations, x)

GI.trait(p::Prepared) = GI.trait(parent(p))
GI.geomtrait(p::Prepared) = GI.geomtrait(parent(p))

GI.isgeometry(::Type{<:Prepared{T}}) where {T} = GI.isgeometry(T)
GI.isfeature(::Type{<:Prepared{T}}) where {T} = GI.isfeature(T)
GI.isfeaturecollection(::Type{<:Prepared{T}}) where {T} = GI.isfeaturecollection(T)

GI.geometry(x::Prepared) = GI.geometry(parent(x))
GI.properties(x::Prepared) = GI.properties(parent(x)) 

for f in (:extent, :crs)
    @eval GI.$f(t::GI.AbstractTrait, x::Prepared) = GI.$f(t, parent(x))
end
for f in (:coordnames, :is3d, :ismeasured, :isempty, :coordinates, :getgeom)
    @eval GI.$f(t::GI.AbstractGeometryTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end

for f in (:x, :y, :z, :m, :coordinates, :getcoord, :ngeom, :getgeom)
    @eval GI.$f(t::GI.AbstractPointTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end
for f in (:npoint, :getpoint, :startpoint, :endpoint, :npoint, :issimple, :isclosed, :isring)
    @eval GI.$f(t::GI.AbstractCurveTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end
for f in (:nring, :getring, :getexterior, :nhole, :gethole, :npoint, :getpoint, :startpoint, :endpoint)
    @eval GI.$f(t::GI.AbstractPolygonTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end
for f in (:npoint, :getpoint, :issimple)
    @eval GI.$f(t::GI.AbstractMultiPointTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
    @eval GI.$f(t::GI.AbstractMultiCurveTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end
for f in (:nring, :getring, :npoint, :getpoint)
    @eval GI.$f(t::GI.AbstractMultiPolygonTrait, geom::Prepared, args...) = GI.$f(t, parent(geom), args...)
end

getpoint(t::GI.AbstractPolyhedralSurfaceTrait, geom::Prepared) = GI.getpoint(t, parent(geom))
isclosed(t::GI.AbstractMultiCurveTrait, geom::Prepared) = GI.isclosed(t, parent(geom))

for f in (:getfeature, :coordinates)
    @eval GI.$f(t::GI.AbstractFeatureTrait, geom::Prepared, args...) = $f(t, parent(geom), args...)
end

# Ambiguity
for T in (:LineTrait, :TriangleTrait, :PentagonTrait, :HexagonTrait, :RectangleTrait, :QuadTrait)
    @eval GI.npoint(t::GI.$T, geom::Prepared) = GI.npoint(t, parent(geom))
end
for T in (:RectangleTrait, :QuadTrait, :PentagonTrait, :HexagonTrait, :TriangleTrait)
    @eval GI.nring(t::GI.$T, geom::Prepared) = GI.nring(t, parent(geom))
end
