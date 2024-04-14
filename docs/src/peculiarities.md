# Peculiarities

## `_True` and `_False` (or `BoolsAsTypes`)

When dispatch can be controlled by the value of a boolean variable, this introduces type instability.  Instead of introducing type instability, we chose to encode our boolean decision variables, like `threaded` and `calc_extent` in `apply`, as types.  This allows the compiler to reason about what will happen, and call the correct compiled method, in a stable way without worrying about 

## What does `apply` return and why?

@rafaqz

## Why do you want me to provide a `target` in set operations?

@skygering

Mainly type stability reasons.
