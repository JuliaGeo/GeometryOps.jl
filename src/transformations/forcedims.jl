function forcexy(geom)
    return GO.apply(GO.GI.PointTrait(), geom) do point
        (GI.x(point), GI.y(point))
    end
end
