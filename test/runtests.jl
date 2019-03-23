using Test
using AlgebraicDependencies

testnum = []
testdep = []

macro num(numbers...)
    push!(testnum, collect(numbers))
end

macro dep(deps...)
    push!(testdep, collect(deps))
end

gsym(n::Int) = [Symbol("v$i") for i in 1:n]

# ------------------------------------------------------------------------------

@num 1/2 2
@dep v1*v2-1

@num 1/2 1 2
@dep v2-1 v1*v3-1

# ------------------------------------------------------------------------------

for (n, d) in zip(testnum, testdep)
    @test dependencies(n; variables=gsym(length(n))) == d
end