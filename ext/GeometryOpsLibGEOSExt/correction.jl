
function (::GO.GEOSCorrection)(::GI.AbstractGeometryTrait, geometry)
    geos_geometry = GI.convert(LG, geometry)
    corrected_geometry = LG.makeValid(geos_geometry)
    @show GI.trait(corrected_geometry)
    if GI.trait(corrected_geometry) isa GI.GeometryCollectionTrait
        corrected_geometry = corrected_geometry[1]
    end
    return corrected_geometry
end