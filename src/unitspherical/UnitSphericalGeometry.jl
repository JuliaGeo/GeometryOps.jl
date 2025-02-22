module UnitSphericalGeometry

import GeoFormatTypes as GFT, GeoInterface as GI

import ..GeometryOps as GO

using LinearAlgebra, ..StaticArrays

"The CRS of all unit spherical geometry.  Note that you have to pass in 3 coordinates to get 3 coordinates!!!"
const UNIT_SPHERICAL_CRS = GFT.WellKnownText(GFT.CRS(), """
GEOCCS["unknown",
    DATUM["unknown",
        SPHEROID["unknown",1,0]
    ],
    PRIMEM["Reference meridian",0],
    UNIT["metre",1,AUTHORITY["EPSG","9001"]],
    AXIS["Geocentric X",OTHER],
    AXIS["Geocentric Y",OTHER],
    AXIS["Geocentric Z",NORTH]]
""")



end