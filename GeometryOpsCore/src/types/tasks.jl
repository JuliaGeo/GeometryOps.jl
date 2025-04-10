
"""
    TaskFunctors(functors, tasks_per_thread)

A struct to hold the functors and tasks_per_thread,
for internal use with `_maptasks` where functions have state
that cannot be shared accross threads, such as `Proj.Transformation`.
"""
struct TaskFunctors{F} 
    functors::F
    tasks_per_thread::Int
end