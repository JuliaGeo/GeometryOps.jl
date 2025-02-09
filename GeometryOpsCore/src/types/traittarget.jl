
#=
## `TraitTarget`

This struct holds a trait parameter or a union of trait parameters.
It's essentially a way to construct unions.
=#

export TraitTarget

"""
    TraitTarget{T}

This struct holds a trait parameter or a union of trait parameters.

It is primarily used for dispatch into methods which select trait levels, 
like `apply`, or as a parameter to `target`.

## Constructors
```julia
TraitTarget(GI.PointTrait())
TraitTarget(GI.LineStringTrait(), GI.LinearRingTrait()) # and other traits as you may like
TraitTarget(TraitTarget(...))
# There are also type based constructors available, but that's not advised.
TraitTarget(GI.PointTrait)
TraitTarget(Union{GI.LineStringTrait, GI.LinearRingTrait})
# etc.
```

"""
struct TraitTarget{T} end
TraitTarget(::Type{T}) where T = TraitTarget{T}()
TraitTarget(::T) where T<:GI.AbstractTrait = TraitTarget{T}()
TraitTarget(::TraitTarget{T}) where T = TraitTarget{T}()
TraitTarget(::Type{<:TraitTarget{T}}) where T = TraitTarget{T}()
TraitTarget(traits::GI.AbstractTrait...) = TraitTarget{Union{map(typeof, traits)...}}()


Base.in(::Trait, ::TraitTarget{Target}) where {Trait <: GI.AbstractTrait, Target} = Trait <: Target
