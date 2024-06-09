using GeometryOps
using Test, SafeTestsets

using GeometryOps.GeoInterface
using GeometryOps.GeometryBasics
using GeoInterface.Extents: Extents
using ArchGDAL
using LibGEOS
using Random, Distributions
using Proj

const GI = GeoInterface
const AG = ArchGDAL
const LG = LibGEOS
const GO = GeometryOps

# Hack LG to have GeometryCollection converts
@eval LG begin
    GI.convert(
        ::Type{GeometryCollection},
        ::GeometryCollectionTrait,
        geom::GeometryCollection;
        context = nothing,
    ) = geom
    function GI.convert(
        ::Type{GeometryCollection},
        ::GeometryCollectionTrait,
        geom;
        context = get_global_context(),
    )
        geometries =
            [
                begin
                    t = GI.trait(g)
                    lg = geointerface_geomtype(t)
                    GI.convert(lg, t, g; context) 
                end
                for g in GI.getgeom(geom)
            ]
        return GeometryCollection(geometries)
    end
end

@testset "GeometryOps.jl" begin
    @safetestset "Primitives" begin include("primitives.jl") end
    # # # Methods
    @safetestset "Angles" begin include("methods/angles.jl") end
    @safetestset "Area" begin include("methods/area.jl") end
    @safetestset "Barycentric coordinate operations" begin include("methods/barycentric.jl") end
    @safetestset "Orientation" begin include("methods/orientation.jl") end
    @safetestset "Centroid" begin include("methods/centroid.jl") end
    @safetestset "DE-9IM Geom Relations" begin include("methods/geom_relations.jl") end
    @safetestset "Distance" begin include("methods/distance.jl") end
    @safetestset "Equals" begin include("methods/equals.jl") end
    # # # Clipping
    @safetestset "Coverage" begin include("methods/clipping/coverage.jl") end
    @safetestset "Cut" begin include("methods/clipping/cut.jl") end
    @safetestset "Polygon Clipping" begin include("methods/clipping/polygon_clipping.jl") end
    # # Transformations
    @safetestset "Embed Extent" begin include("transformations/extent.jl") end
    @safetestset "Reproject" begin include("transformations/reproject.jl") end
    @safetestset "Flip" begin include("transformations/flip.jl") end
    @safetestset "Simplify" begin include("transformations/simplify.jl") end
    @safetestset "Segmentize" begin include("transformations/segmentize.jl") end
    @safetestset "Transform" begin include("transformations/transform.jl") end
    @safetestset "Geometry correction" begin 
        include("transformations/correction/geometry_correction.jl")
        include("transformations/correction/closed_ring.jl") 
        include("transformations/correction/intersecting_polygons.jl")
    end
    # Extensions
    @safetestset "FlexiJoins" begin include("extensions/flexijoins.jl") end
end
