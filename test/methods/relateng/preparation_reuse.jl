using Test, Random
import GeometryOps as GO, GeoInterface as GI

@testset "Spherical preparation reuses transient bounds conservatively" begin
    function encloses(a,b)
        all(getproperty(a,d)[1] <= getproperty(b,d)[1] + 2e-14 &&
            getproperty(a,d)[2] >= getproperty(b,d)[2] - 2e-14 for d in (:X,:Y,:Z))
    end
    rng = MersenneTwister(59)
    shells = [
        [(0.,80.), (120.,80.), (240.,80.), (0.,80.)],
        [(0.,0.), (30.,0.), (15.,30.), (0.,0.)],
        [(170.,-10.), (-170.,-10.), (-170.,10.), (170.,10.), (170.,-10.)],
        [(2.,0.,0.), (0.,3.,0.), (0.,0.,4.), (2.,0.,0.)],
    ]
    append!(shells, [[(lon+5cos(t),lat+4sin(t)) for t in range(0,2pi;length=25)]
        for (lon,lat) in [(rand(rng)*320-160,rand(rng)*120-60) for _ in 1:20]])
    for coords in shells, reversed in (false,true), closed in (false,true), oriented in (false,true)
        c = reversed ? reverse(coords) : coords
        c = closed ? c : c[1:end-1]
        shell = GI.LinearRing(c)
        poly = GI.Polygon([shell])
        m = GO.Spherical(; oriented)
        original = GO.rk_interaction_bounds(m,poly)
        cached = GO._relate_cache_extents(m,poly)
        @test encloses(cached.extent, original)
        normalized = GO._ring_usp(shell)
        @test all(eachindex(normalized)) do i
            a, b = normalized[i], normalized[mod1(i+1,length(normalized))]
            all((0.,0.25,0.5,0.75,1.)) do t
                q = GO.UnitSpherical.slerp(a,b,t)
                all(getproperty(cached.extent,d)[1] <= q[j] <= getproperty(cached.extent,d)[2]
                    for (j,d) in enumerate((:X,:Y,:Z)))
            end
        end
        @test GO._relate_cache_extents(m,cached) === cached
        @test GI.coordinates(GI.getexterior(cached)) == GI.coordinates(shell)
    end
    shell = GI.LinearRing([(0.,40.), (20.,40.), (20.,55.), (0.,55.), (0.,40.)])
    hole = GI.LinearRing([(8.,45.), (12.,45.), (12.,49.), (8.,49.), (8.,45.)])
    stray = GI.LinearRing([(80.,0.), (85.,0.), (85.,5.), (80.,5.), (80.,0.)])
    for rings in ([shell,hole], [shell,stray]), oriented in (false,true)
        m = GO.Spherical(;oriented)
        poly = GI.Polygon(rings)
        cached = GO._relate_cache_extents(m,poly)
        @test encloses(cached.extent,GO.rk_interaction_bounds(m,poly))
        for r in GI.getring(cached)
            @test encloses(cached.extent,r.extent)
        end
        # Pre-stamped rings retain the original fallback semantics.
        cached_rings = GI.Polygon([GO._relate_cache_extents(m,r) for r in rings])
        @test encloses(GO._relate_cache_extents(m,cached_rings).extent,GO.rk_interaction_bounds(m,cached_rings))
    end
    m = GO.Spherical()
    stamped_hole_poly = GI.Polygon([shell, GO._relate_cache_extents(m, hole)])
    @test GO._relate_cache_extents(m, stamped_hole_poly).extent == GO.rk_interaction_bounds(m, stamped_hole_poly)
    @test_throws ArgumentError GO._relate_cache_extents(m,GI.LineString([(1.,0.,0.),(-1.,0.,0.)]))
    @test_throws ArgumentError GO._relate_cache_extents(m,GI.Polygon([GI.LinearRing([(1.,0.,0.),(-1.,0.,0.),(0.,1.,0.),(1.,0.,0.)])]))
    exactring = GI.LinearRing([(2.,0.,0.),(0.,3.,0.),(0.,0.,4.),(2.,0.,0.)])
    cachedring = GO._relate_cache_extents(m, exactring)
    @test Tuple(GO._ring_kernel_pts(cachedring)[1]) == (2.,0.,0.)
    @test Tuple(GO._ring_usp(cachedring)[1]) == (1.,0.,0.)
end
