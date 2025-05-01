# How to write fast code

In the GeoInterface ecosystem and GeometryOps specifically, there are a few tricks that can help you keep your code fast and allocation free.

## Always propagate compile-time information

The first time you call `trait` should be the last time you call `trait` on that geometry.  
Otherwise - propagate that trait down the stack!

If you don't, then the compiler loses track of it, and when it finds it again, it has to allocate to perform a dynamic dispatch.  
This is pretty slow and can cause a 3x (or much much larger) slowdown in your code.

Things like the [`Applicator`](@ref)s and especially the [`ApplyWithTrait`](@ref) applicator can help here.

Similarly, you'll notice a pattern where we pass a floating point type down the chain.  This is done for type stability as well.
If GeoInterface gets a `coordtype` in future then it'll default to `float(coordtype(geom))`, but for now we fix it at f64 and let
the user change it if they want.  This lets us avoid all the issues with "oh but I have a float32 geometry or a bigfloat 
geometry or something".

## Try not to allocate unless necessary

There are a lot of algorithms that seem simple to implement with some `collect`s.  Try to skip that if possible, and use
GeoInterface constructs like `getgeom`, `getpoint`, and `getring`, which are faster anyway.

## Analyse your code using Julia tools

[**ProfileView.jl**](https://github.com/timholy/ProfileView.jl) and [**Cthulhu.jl**](https://github.com/Cthulhu.jl) work together very well to diagnose
and fix type instability.  [**JET.jl**](https://github.com/aviatesk/JET.jl) is also good here.

[**TimerOutputs.jl**](https://github.com/KristofferC/TimerOutputs.jl) is excellent for characterizing where your time is being spent 
and which parts of your function you should focus on optimizing.  Always use TimerOutputs before hyperoptimizing - you don't usually 
want to halve the cost of a function which contributes 1% of your runtime!

## Use statically sized, immutable types where you can

Static, immutable types are very good because they can be inlined and do not allocate.  
But this isn't a taboo against mutables by any means.  Sometimes rolling your own stack (which allocates) is substantially faster than
recursion (which technically doesn't).

If you have a type which you don't know the size of, and which you believe is completely random and unpredictable at compile time,
pay the cost and make it a vector instead of forcing type instability.  This applies to tuples etc.

