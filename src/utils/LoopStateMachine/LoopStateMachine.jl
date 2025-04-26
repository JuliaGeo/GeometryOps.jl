"""
    LoopStateMachine

Utilities for returning state from functions that run inside a loop.

This is used in e.g clipping, where we may need to break or transition states.

The main entry point is to return an [`Action`](@ref) from a function that 
is wrapped in a `@controlflow f(...)` macro in a loop.  When a known `Action`
(currently, `:continue`, `:break`, `:return`, or `:full_return` actions) is returned, 
it is processed by the `@controlflow` macro, which allows the function to break out of the loop
early, continue to the next iteration, or return a value, basically a way to provoke syntactic
behaviour from a function called from a inside a loop, where you do not have access to that loop.

## Example

```julia
```
"""
module LoopStateMachine

export Action, @controlflow

import ..GeometryOps as GO

const ALL_ACTION_DESCRIPTIONS = """
- `:continue`: continue to the next iteration of the loop.  
  This is the `continue` keyword in Julia.  The contents of the action are not used.
- `:break`: break out of the loop.  
  This is the `break` keyword in Julia.  The contents of the action are not used.
- `:return`: cause the function executing the loop to return with the wrapped value.
- `:full_return`: cause the function executing the loop to return `Action(:full_return, x)`.  
  This is very useful to terminate recursive funtions, like tree queries terminating after you 
  have found a single intersecting segment.
"""

"""
    Action(name::Symbol, [x])

Create an `Action` with the name `name` and optional contents `x`.

`Action`s are returned from functions wrapped in a `@controlflow` macro, which
does something based on the return value of that function if it is an `Action`.

## Available actions

$ALL_ACTION_DESCRIPTIONS
"""
struct Action{T}
    name::Symbol
    x::T
end

Action() = Action{Nothing}(:unnamed, nothing)
Action(x::T) where T = Action{T}(:unnamed, x)
Action(x::Symbol) = Action(x, nothing)

function Base.show(io::IO, action::Action{T}) where T
    print(io, "Action")
    print(io, "(:$(action.name)")
    if isnothing(action.x)
        print(io, ")")
    else
        print(io, ", ",action.x, ")")
    end
end

struct UnrecognizedActionException <: Base.Exception
    name::Symbol
end

function Base.showerror(io::IO, e::UnrecognizedActionException)
    print(io, "Unrecognized action: ")
    printstyled(io, e.name; color = :red, bold = true)
    println(io, ".")
    println(io, "Valid actions are:")
    println(io, ALL_ACTION_DESCRIPTIONS)
end

# We exclude the macro definition from code coverage computations,
# because I know it's tested but Codecov doesn't seem to think so.
# COV_EXCL_START
"""
    @controlflow f(...)

Process the result of `f(...)` and return the result if it's not an `Action`(@ref LoopStateMachine.Action).
    
If it is an `Action`, then process it according to the following rules, and throw an error if it's not recognized.
`:continue`, `:break`, `:return`, or `:full_return` are valid actions.

$ALL_ACTION_DESCRIPTIONS

!!! warning
    Only use this inside a loop, otherwise you'll get a syntax error, especially if you use `:continue` or `:break`.

## Examples
"""
macro controlflow(expr)
    varname = gensym("loop-state-machine-returned-value")
    return quote
        $varname = $(esc(expr))
        if $varname isa Action
            if $varname.name == :continue
                continue
            elseif $varname.name == :break
                break
            elseif $varname.name == :return
                return $varname.x
            elseif $varname.name == :full_return
                return $varname
            else
                throw(UnrecognizedActionException($varname.name))
            end
        else
            $varname
        end
    end
end
# COV_EXCL_STOP

# You can define more actions as you desire.

end