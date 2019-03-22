module AlgebraicDependencies

export dependencies

using SymEngine
using SymPy
using Polynomials
using Hecke
using Nemo
using LinearAlgebra
using ContinuedFractions
using Singular

function dependencies(roots::Vector{Basic}; variables=Basic[])
    # println("Computing dependencies between ", roots)
    if length(roots) < 2
        return nothing
    end
    lattice = findrelations(roots)
    # if isempty(variables)
    #     variables = symset("v", length(lattice))
    # elseif length(variables) != ncols(lattice)
    #     throw("Number of variables does not match number of columms. Got $(length(variables)), need $(ncols(lattice))")
    # end
    return ideal(lattice, variables)
end

is_rational(x::Basic) = SymEngine.BasicType(x) isa SymEngine.BasicType{Val{:Integer}} || SymEngine.BasicType(x) isa SymEngine.BasicType{Val{:Rational}}

function minpoly(r::Basic, z::Symbol)
    s = Sym(z)
    p = SymPy.minpoly(Sym(string(r)), s)
    cs = SymPy.coeffs(p, s)
    bs = map(Basic ∘ string, cs)
    Basic(string(p)), Polynomials.Poly(bs)
end

function common_number_field(roots::Vector{Basic})
    v = :a
    minpolys = Basic[]
    ps = Polynomials.Poly[]
    for r in roots
        mp, p = minpoly(r, v)
        if is_rational(r)
            push!(minpolys, r)
        else
            push!(minpolys, mp)
        end
        push!(ps, p)
    end
    indices = findall(is_rational, minpolys)
    ratroots = [minpolys[i] for i in indices] 
    deleteat!(minpolys, indices)
    mps = [convert(Expr, b) for b in minpolys]
    if isempty(mps)
        return 1, roots, ps 
    end

    deleteat!(roots, indices)
    R, _ = PolynomialRing(QQ, string(v))
    K, g = number_field(map(R, mps))
    @debug "Number field" K
    S, mS = simple_extension(K)
    @debug "Simple extension" S
    gs = [SymEngine.lambdify(Basic(replace(string(mS\(x)), "//"=>"/"))) for x in g]
    rs = [f isa Function ? f(r) : Basic(f) for (f, r) in zip(gs, roots)]
    for (i, ind) in enumerate(indices)
        splice!(rs, ind:ind-1, ratroots[i])
    end
    Hecke.degree(S), rs, ps
end

function masser_bound(degree::Int, roots::Vector{Basic}, mpolys::Vector{Polynomials.Poly})
    k = length(roots)
    # assume all roots belong to the same number field
    d = degree
    # h = maximum height of the a[i]. The height of an algebraic number is the 
    # sum of the degree and the binary length of all coefficients in the 
    # defining equation over Q
    h = 0
    for p in mpolys
        cs = filter(x->!iszero(x), Polynomials.coeffs(p))
        h0 = ceil(Polynomials.degree(p) + sum([log(abs(c)) for c in cs]))
        if h0 > h
            h = h0
        end
    end
    @debug "Parameters for Masser bound" d h k
    return ceil(d^2 * (4*h*k*d* (log(2+d)/log(log(2+d)))^3)^(k-1) + 1)
end

function ideal(m::Matrix{BigInt}, x::Array{Basic,1})
    if isempty(m)
        return []
    end

    y = [Basic("v$i") for i in 1:length(x)]
    base = Basic[]
    for i in 1:nrows(m)
        r = 1
        for j in 1:ncols(m)
            exp = m[i,j]
            v = exp > 0 ? x[j] : y[j]
            r *= v^abs(exp)
        end
        base = [base; r - 1]
    end

    inv = [xi*yi - 1 for (xi,yi) in zip(x,y)]
    base = [base; inv]
    
    R, g = Singular.PolynomialRing(Singular.QQ, map(string, [y; x]))
    sbasis = [R(convert(Expr, b)) for b in base]
    I = Singular.Ideal(R, sbasis)
    return Singular.eliminate(I, g[1:length(y)]...)
end


function findrelations(roots::Vector{Basic})
    # first treat zeros in the root list
    # zeros = find(x -> x == 0, roots)
    # if length(zeros) == length(roots)
    #     return []
    # end
    # if !isempty(zeros)
    #     B = findrelations(filter(x->x!=0, roots))
    #     # TODO: insert new dimensions
    #     return B
    # end

    # TODO: common number field does not work as expected
    degree, roots, mpolys = common_number_field(roots) # TODO: does nothing if the roots belong to the same field
    @debug "Common number field" degree roots mpolys
    rs = [convert(Complex{BigFloat}, x) for x in roots]
    @debug "Complex roots" rs

    relog = [real(log(x)) for x in rs]
    imlog = [imag(log(x)) for x in rs]
    imlog = [imlog; 2*BigFloat(pi)]

    # Replace implicit zeros by explicit ones
    for i in 1:length(rs)
        z = rs[i]

        # abs(z) == 1
        if abs(abs(z) - 1) < 0.1 && abs(z) == 1
            relog[i] = 0
        end

        # z is real and z >= 0
        if isreal(z) && isreal(sqrt(z))
            imlog[i] = 0
        else
            # TODO: try harder: If[ Element[RootReduce[z], Reals] && Element[Sqrt[RootReduce[z]], Reals], imLog[[i]] = 0 ];
        end
    end
    @debug "reLog and imLog" relog imlog


    # Compute a bound for the coefficients of the generators
    bound = Int(masser_bound(degree, roots, mpolys))
    @debug "Masser bound" bound

    m = eye(Int, length(roots)) # identity matrix

    # Successively refine the approximation until only valid generators are returned
    level = Int(ceil(log2(bound) + 1))
    while prod(Bool[check_relation(roots, m[i,:]) for i in 1:size(m)[1]]) == 0
        m1 = getbasis(relog, level, bound)
        m2 = getbasis(imlog, level, bound)
        m = z_module_intersect(m1, m2[:,1:end-1])
        level = level + 1
    end

    return m
end

hnf_with_transform(m::Matrix{Int}) = Matrix{Int}.(Nemo.hnf_with_trafo(matrix(FlintZZ, m)))

function lattice_divide(l::Matrix{Int}, d::Int)
    n = ncols(l)
    return z_module_intersect(l, d*eye(Int, n)) / d
end

function check_relation(a::Vector{Basic}, e::Vector{<:Integer})
    return SymEngine.expand(prod([ax^ex for (ax, ex) in zip(a,e)])) == 1
end

function convergent(x, n)
    cf = ContinuedFraction(x)
    co = convergents(cf)
    return iterate(co, n)[1]
end

nrows(a::Matrix{<:Any}) = size(a)[1]
ncols(a::Matrix{<:Any}) = size(a)[2]

function getbasis(l::Vector{BigFloat}, level::Int, bound::Int)
    n = length(l)

    # First treat zero elements in l as special case
    zpos = findall(iszero, l)
    if length(zpos) == length(l)
        return eye(Int, n)
    end

    if length(zpos) > 0
        ll = deleteat!(copy(l), zpos)
        b = getbasis(ll, level, bound) # basis for nonzero numbers
        zvec = zeros(nrows(b), 1)
        # Insert new dimensions
        for pos in zpos
            b = hcat(b[:,1:pos-1], zvec, b[:,pos:end])
        end
        # Add unit vectors at the zero positions
        b = vcat(b, eye(Int, n)[zpos,:])
        return Matrix{Int}(b)
    end

    # Now we assume that l is a list of nonzero real numbers
    c1 = [convergent(x, level) for x in l]
    c2 = [convergent(x, level+1) for x in l]

    e = [1//denominator(x1*x2) for (x1, x2) in zip(c1, c2)]
    # Refine the approximation such that all errors are smaller than the smallest
    # element of l in absolute value
    lev = level + 1
    while length(filter(x -> maximum(e) >= abs(x), c1)) > 0 && lev < level + 5
        ex = findall(x->x==maximum(e), e) # indices with greates error

        lev++
        for i in 1:length(ex)
            j = ex[i]
            c1[j] = c2[j]
            c2[j] = convergent(l[j], lev)
            e[j] = c1[j] == l[j] ? 0 : 1/denominator(c1[j]*c2[j])
        end
    end

    # now: max e[i] < min |c1[i]|
    # this bound guarantees that generators whose norm is at most bound will
    # appear in the LLL-reduced basis
    minc1 = minimum([abs(c) for c in c1])
    maxe = maximum([abs(c) for c in e])
    d = BigInt(ceil(2^((length(c1) - 1)/2)*bound/(minc1 - maxe)))
    identity = eye(Int, n)
    row = c1 * d
    b = hcat(identity, row)
    b = lll(b)
    # Vectors whose right hand side is greater than the errors permit can be 
    #   discarded; they cannot correspond to integer relations.
    # b = vcat([b[i,:] for i in 1:nrows(b) ])
    # TODO: find better way to filter rows
    res = Matrix{Rational{Int}}(undef, 0,n+1)
    for i in 1:nrows(b)
        if (abs(b[i,:][end]) <= d*abs(dot(b[i,1:n],e)))
            res = vcat(res, transpose(b[i,:]))
        end
    end

    # all remaining vectors are returned as candidates
    # TODO: result should be integer matrix?
    return Matrix{Int}(res[:,1:end-1])
end

function z_nullspace(matrix::Matrix{Int})
    h, t = hnf_with_transform(matrix)
    t = t * -1
    # println("HNF: $(h) | $(t)")

    # kernel is generated by the rows of t that correspond to zero rows in h
    zvec = zeros(size(h, 2))

    # TODO: find better way to filter zero vectors
    res = Matrix{Int}(undef, 0, ncols(t))
    for i in 1:nrows(t)
        if iszero(h[i,:])
            res = vcat(res, transpose(t[i,:]))
        end
    end
    return res
    # return [t[i,:] for i in 1:size(h, 1) if h[i,:] == zvec]
end

function z_module_intersect(base1::Matrix{Int}, base2::Matrix{Int})
    if isempty(base1) || isempty(base2)
        return []
    end

    sol = z_nullspace(vcat(base1, base2))

    if isempty(sol)
        return []
    end

    m1 = transpose(base1)
    m2 = transpose(sol[:, 1:nrows(base1)])
    return lll(collect(transpose(m1 * m2)))
end

function clear_denom(a::Matrix{Rational{BigInt}})
    d = lcm(denominator.(a))
    return a*d, d
end

import Nemo.lll

function lll(a::Matrix{Rational{BigInt}})
    m, d = clear_denom(a)
    m = numerator.(m)
    m = Matrix{BigInt}(Nemo.lll(matrix(FlintZZ, m)))
    return m // d
end

lll(m::Matrix{Int}) = Matrix{BigInt}(Nemo.lll(matrix(FlintZZ, m)))

eye(::Type{T}, n) where {T} = Matrix{T}(I, n, n)

function (R::FmpqPolyRing)(p::Expr)
    v = gen(R)
    vs = [:($(Symbol(string(v))) = $v)]
    q = quote
        let $(vs...)
            $(p)
        end
    end
    eval(q)
end

function (R::Singular.PolyRing)(p::Expr)
    vs = gens(R)
    vs = [:($(Symbol(string(v))) = $v) for v in vs]
    q = quote
        let $(vs...)
            $(p)
        end
    end
    eval(q)
end

function (K::NfAbsNS)(a::fmpq_poly)
    q, w = divrem(a, K.pol)
    z = NfAbsNSElem(w)
    z.parent = K
    return z
end

end # module
