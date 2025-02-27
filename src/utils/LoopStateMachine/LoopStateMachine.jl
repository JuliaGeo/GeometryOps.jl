"""
    LoopStateMachine

Utilities for returning state from functions that run inside a loop.

This is used in e.g clipping, where we may need to break or transition states.

The main entry point is to return an [`Action`](@ref) from a function that 
is wrapped in a `@controlflow f(...)` macro in a loop.  When a known `Action`
(currently, `Continue`, `Break`, or `Return`) is returned, it is processed
by the `@controlflow` macro, which allows the function to break out of the loop
early, continue to the next iteration, or return a value, without being 
syntactically inside the loop.
"""
module LoopStateMachine

struct Action{name, T}
    x::T
end

Action{name}() where name = Action{name, Nothing}(nothing)
Action{name}(x::T) where name where T = new{name, T}(x)
Action(x::T) where T = Action{:unnamed, T}(x)

Action{name, Nothing}() where name = Action{name, Nothing}(nothing)

function Base.show(io::IO, action::Action{name, T}) where {name, T}
    print(io, "Action ", name)
    if isnothing(action.x)
        print(io, "()")
    else
        print(io, "(", action.x, ")")
    end
end


# Some common actions
"""
    Break()

Break out of the loop.
"""
const Break = Action{:Break, Nothing}

"""
    Continue()

Continue to the next iteration of the loop.
"""
const Continue = Action{:Continue, Nothing}

"""
    Return(x)

Cause the function executing the loop to return.  Use with great caution!
"""
const Return = Action{:Return}

"""
    @controlflow f(...)

Process the result of `f(...)` and return the result if it's not a [`Continue`](@ref), [`Break`](@ref), or [`Return`](@ref) [`Action`](@ref).

- `Continue`: continue to the next iteration of the loop.
- `Break`: break out of the loop.
- `Return`: cause the function executing the loop to return with the wrapped value.

!!! warning
    Only use this inside a loop, otherwise you'll get a syntax error!
"""
macro controlflow(expr)
    varname = gensym("lsm-f-ret")
    return quote
        $varname = $(esc(expr))
        if $varname isa Continue
            continue
        elseif $varname isa Break
            break
        elseif $varname isa Return
            return $varname.x
        else
            $varname
        end
    end
end

# You can define more actions as you desire.

end