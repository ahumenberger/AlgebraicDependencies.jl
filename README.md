# AlgebraicDependencies.jl

A Julia package for computing algebraic dependencies among exponential sequences.

```julia
julia> @deps 2 1/2
1-element Array{Expr,1}:
 :(v1 * v2 - 1)

julia> @deps (1-sqrt(5))/2 (1+sqrt(5))/2
1-element Array{Expr,1}:
 :(v1 ^ 2 * v2 ^ 2 - 1)
```
