using Test
using GeometryOps.LoopStateMachine: @controlflow, Action

@testset "Continue action" begin
    count = 0
    f(i) = begin
        count += 1
        if i == 3
            return Action(:continue)
        end
        count += 1
    end
    for i in 1:5
        @controlflow f(i)
    end
    @test count == 9 # Adds 1 for each iteration, but skips second +1 on i=3
end

@testset "Break action" begin
    count = 0
    function f(i)
        count += 1
        if i == 3
            return Action(:break)
        end
        count += 1
    end
    for i in 1:5
        @controlflow f(i)
    end
    @test count == 5 # Counts up to i=3, adding 2 for i=1,2 and 1 for i=3
end

@testset "Return action" begin
    f(i) = for j in 1:3
        i == j && @controlflow Action(:return, i)
    end 
    @test f(1) == 1
    @test f(2) == 2
    @test f(3) == 3
end

@testset "Full return action" begin
    f(i) = for j in 1:3
        i == j && @controlflow Action(:full_return, i)
    end
    @test f(1) == Action(:full_return, 1)
    @test f(2) == Action(:full_return, 2)
    @test f(3) == Action(:full_return, 3)
end

@testset "Return value" begin
    results = Int[]
    for i in 1:3
        val = @controlflow begin
            i * 2
        end
        push!(results, val)
    end
    @test results == [2, 4, 6]
end

@testset "Show" begin
    @test sprint(print, Action(:continue)) == "Action(:continue)"
    @test sprint(print, Action(:break)) == "Action(:break)"
    @test sprint(print, Action(:return, 1)) == "Action(:return, 1)"
    @test sprint(print, Action(:full_return, 1)) == "Action(:full_return, 1)"
end

@testset "Unnamed action" begin
    @test sprint(print, Action()) == "Action(:unnamed)"
    @test sprint(print, Action(1)) == "Action(:unnamed, 1)"
    @test sprint(print, Action(:x)) == "Action(:x)"
end

