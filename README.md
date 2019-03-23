# AlgebraicDependencies.jl

A Julia package for computing algebraic dependencies among exponential sequences. Let `r_1,...,r_s` be algebraic numbers. We can compute a basis for the ideal `I` containing all algebraic dependencies among `r_1^n,...,r_s^n`, i.e. for every `p` in `I` we have `p(r_1^n,...,r_s^n) = 0` for all natural numbers `n`.

```julia
julia> @deps 2 1/2
1-element Array{Expr,1}:
 :(v1 * v2 - 1)

julia> @deps (1-sqrt(5))/2 (1+sqrt(5))/2
1-element Array{Expr,1}:
 :(v1 ^ 2 * v2 ^ 2 - 1)
```
