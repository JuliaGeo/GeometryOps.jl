import GeoInterface as GI, GeometryOps as GO, LibGEOS as LG
using SortTileRecursiveTree, Tables, DataFrames
using NaturalEarth # data source

using GeoMakie # enable plotting
# We have to monkeypatch GeoMakie.to_multipoly for geometry collections:
function GeoMakie.to_multipoly(::GI.GeometryCollectionTrait, geom)
    geoms = collect(GI.getgeom(geom))
    poly_and_multipoly_s = filter(x -> GI.trait(x) isa GI.PolygonTrait || GI.trait(x) isa GI.MultiPolygonTrait, geoms)
    if isempty(poly_and_multipoly_s)
        return GeometryBasics.MultiPolygon([GeometryBasics.Polygon(Point{2 + GI.hasz(geom) + GI.hasm(geom), Float64}[])])
    else
        final_multipoly = reduce((x, y) -> GO.union(x, y; target = GI.MultiPolygonTrait()), poly_and_multipoly_s)
        return GeoMakie.to_multipoly(final_multipoly)
    end
end


function _geom_vector(object)
    if GI.trait(object) isa GI.FeatureCollectionTrait
        return GI.geometry.(GI.getfeature(object))
    elseif Tables.istable(object)
        return Tables.getcolumn(object, first(GI.geometrycolumns(object)))
    elseif object isa AbstractVector
        return object
    else
        error("Given object was neither a FeatureCollection, Table, nor AbstractVector.  Could not extract a vector of geometries.  Type was $(typeof(object))")
    end
end

function adjacency_matrix(weight_f, predicate_f, source, target; self_intersection = false, use_strtree = false)
    @info "Config: use_strtree = $use_strtree, self_intersection = $self_intersection"

    source_geoms = _geom_vector(source)
    target_geoms = _geom_vector(target)

    local target_rtree
    if use_strtree
        # Create an STRTree on the target geometries
        target_rtree = SortTileRecursiveTree.STRtree(target_geoms)
    end

    # create an adjacency matrix
    adjacency_matrix = zeros(length(source_geoms), length(target_geoms))

    if use_strtree
        # loop over the source tiles and target tiles and fill in the adjacency matrix
        # only call weight_f on those geometries which pass:
        # (a) the STRTree test and 
        # (b) the predicate_f test
        # WARNING: sometimes, STRTree may not provide correct results.
        # Maybe use LibSpatialIndex instead?
        for (i, source_geom) in enumerate(source_geoms)
            targets_in_neighbourhood = SortTileRecursiveTree.query(target_rtree, source_geom)
            for target_index in targets_in_neighbourhood
                if predicate_f(source_geom, target_geoms[target_index])
                    adjacency_matrix[i, target_index] = weight_f(source_geom, target_geoms[target_index])
                end
            end
        end
    else
        # loop over the source tiles and target tiles and fill in the adjacency matrix
        # only call weight_f on those geometries which pass:
        # (a) the predicate_f test
        # TODO: potential optimizations:
        # - if weight_f is symmetric, then skip the computation if that index of the matrix
        #   is already nonzero, i.e., full
        for (i, source_geom) in enumerate(source_geoms)
            for (j, target_geom) in enumerate(target_geoms)
                if !self_intersection && i == j # self intersection
                    continue
                end
                if predicate_f(source_geom, target_geom)
                    adjacency_matrix[i, j] = weight_f(source_geom, target_geom)
                end
            end
        end
    end

    if self_intersection
        # fill in the identity line
        for i in 1:length(source_geoms)
            adjacency_matrix[i, i] = 1.0
        end
    end

    return adjacency_matrix
end

# basic test using NaturalEarth US states

# get the US states as a GeoInterface FeatureCollection
us_states = DataFrame(naturalearth("admin_1_states_provinces", 110)) # 110m is only US states but 50m is all states, be careful then about filtering.
# We also have to make all geometries valid so that they don't cause problems for the predicates!
us_states.geometry = LG.makeValid.(GI.convert.((LG,), us_states.geometry))

@time adjmat = adjacency_matrix(LG.touches, us_states, us_states) do source_geom, target_geom
    return 1.0
    # this could be some measure of arclength, for example.
end

f, a, p = poly(us_states.geometry |> GeoMakie.to_multipoly; color = fill(:black, size(us_states, 1)))

record(f, "adjacency_matrix.mp4", 1:size(us_states, 1); framerate = 1) do i
    a.title = "$(us_states[i, :name])"
    colors = fill(:gray, size(us_states, 1))
    colors[(!iszero).(view(adjmat, :, i))] .= :yellow
    colors[i] = :red
    p.color = colors
end
