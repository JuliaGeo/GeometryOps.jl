#=
# TaskFunctors

Many functions a user might run through `apply` call into C.  But what do you do
when the C function you are calling is not threadsafe, or uses a reentrant/context
based API?

The answer is to wrap the function in TaskFunctors, which will use an approximation of
task local storage to ensure that each task calls its own copy of the function or C object
you are invoking.

The primary application of this is Proj `reproject`, and Proj has a reentrant API based on context
objects.  So, to make this work, we need to create a new context object per task.  We do this and 
pass the vector of Proj transformation objects to `TaskFunctors`, which `apply` and `_maptasks`
dispatch on to behave correctly in this circumstance.

```@docs; canonical=false
TaskFunctors
```
=#
"""
    TaskFunctors(functors)

A struct to hold the functors,
for internal use with `_maptasks` where functions have state
that cannot be shared accross threads, such as `Proj.Transformation`.

`functors` must be an array or tuple of functors, one per task.

The number of tasks is the number of elements in `functors`.
"""
struct TaskFunctors{F} 
    functors::F
end