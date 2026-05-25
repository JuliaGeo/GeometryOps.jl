# # JTS NG algorithm types

export RelateNG, OverlayNG

abstract type NGPlanarAlgorithm <: GeometryOpsCore.Algorithm{Planar} end

best_manifold(::NGPlanarAlgorithm, input) = Planar()
rebuild(alg::T, ::Planar) where {T <: NGPlanarAlgorithm} = alg
rebuild(alg::T, ::AutoManifold) where {T <: NGPlanarAlgorithm} = alg
function rebuild(alg::T, ::M) where {T <: NGPlanarAlgorithm, M <: Manifold}
    throw(
        GeometryOpsCore.WrongManifoldException{M, Planar, T}(
            "The JTS NG ports are planar algorithms.",
        ),
    )
end

"""
    RelateNG(; boundary_node_rule = Mod2BoundaryNodeRule(), prepared = false)

Algorithm marker for the JTS RelateNG port.  Predicate methods will dispatch on
this type as the engine is ported.
"""
struct RelateNG{R <: BoundaryNodeRule} <: NGPlanarAlgorithm
    boundary_node_rule::R
    prepared::Bool
end

RelateNG(; boundary_node_rule::BoundaryNodeRule = Mod2BoundaryNodeRule(), prepared::Bool = false) =
    RelateNG(boundary_node_rule, prepared)
RelateNG(::Planar; kwargs...) = RelateNG(; kwargs...)

"""
    NoPrecisionModel()

Placeholder precision model for full-precision OverlayNG operation.
"""
struct NoPrecisionModel end

"""
    OverlayNG(; strict = false, area_result_only = false, optimized = true,
              precision_model = NoPrecisionModel())

Algorithm marker and configuration holder for the JTS OverlayNG port.
Constructive overlay methods will dispatch on this type as the engine is ported.
"""
struct OverlayNG{P} <: NGPlanarAlgorithm
    strict::Bool
    area_result_only::Bool
    optimized::Bool
    precision_model::P
end

function OverlayNG(;
    strict::Bool = false,
    area_result_only::Bool = false,
    optimized::Bool = true,
    precision_model = NoPrecisionModel(),
)
    return OverlayNG(strict, area_result_only, optimized, precision_model)
end
OverlayNG(::Planar; kwargs...) = OverlayNG(; kwargs...)
