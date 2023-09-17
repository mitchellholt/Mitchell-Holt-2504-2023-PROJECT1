#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity)
such that their product of the list (mod p) is f. Irreducibles are fixed points
on the function factor.
"""
function factor(f :: PolynomialModP) :: Vector{Tuple{PolynomialModP, Int}}
    #Cantor Zassenhaus factorization
    degree(f) <= 1 && return [(f, 1)]

    # make f square-free
    squares_poly = gcd(f, derivative(ff))
    # @show squares_poly, typeof(squares_poly)
    ff = f ÷ squares_poly
    # @show "after square free:", ff, typeof(ff)

    # make f monic
    old_coeff = leading(ff).coeff
    ff = ff ÷ old_coeff
    # @show "after monic:", ff, typeof(ff)

    dds = dd_factor(ff)

    ret_val = Tuple{PolynomialModP, Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k)
        sp = map(p -> p ÷ leading(p).coeff, sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f, mp)))
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f).coeff * one(f), 1))

    return ret_val
end

function factor(f :: P, prime :: Int) :: Vector{Tuple{P, Int}} where P <: Polynomial
    #Cantor Zassenhaus factorization

    f_modp = mod(f, prime)
    degree(f_modp) <= 1 && return [(f_modp, 1)]

    # make f primitive
    ff = prim_part(f_modp)(prime)
    # @show "after prim:", ff, typeof(ff)

     # make f square-free
    squares_poly = gcd(f, derivative(ff), prime)
    # @show squares_poly, typeof(squares_poly)
    ff = (ff ÷ squares_poly)(prime)
    # @show "after square free:", ff, typeof(ff)

    # make f monic
    old_coeff = leading(ff).coeff
    ff = (ff ÷ old_coeff)(prime)        
    # @show "after monic:", ff, typeof(ff)

    dds = dd_factor(ff, prime)

    ret_val = Tuple{P ,Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k, prime)
        sp = map((p)->(p ÷ leading(p).coeff)(prime),sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, multiplicity(f_modp,mp,prime)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (leading(f_modp).coeff* one(P), 1) )

    return ret_val
end

"""
Expand a factorization.
"""
function expand_factorization(factorization :: Vector{Tuple{P, Int}}) where P <: Polynomial
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f
"""
function multiplicity(f :: PolynomialModP, g :: PolynomialModP) :: Int
    degree(gcd(f, g)) == 0 && return 0
    return 1 + multiplicity(f ÷ g, g)
end

function multiplicity(f :: P, g :: P, prime :: Int) :: Int where P <: Polynomial
    degree(gcd(f, g, prime)) == 0 && return 0
    return 1 + multiplicity((f ÷ g)(prime), g, prime)
end


"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a
product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2,
such that the product of the list (mod `prime`) is equal to `f` (mod `prime`).
"""
function dd_factor(f :: PolynomialModP) :: Array{PolynomialModP}
    x = x_poly(PolynomialModP, f.prime)
    w = deepcopy(x)
    g = Array{PolynomialModP}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem(w^prime, f)
        g[k] = gcd(w - x, f)
        f = f ÷ g[k]
    end

    #edge case for final factor
    f != one(f) && push!(g,f)

    return g
end

function dd_factor(f :: P, prime :: Int) :: Array{P} where P <: Polynomial
    x = x_poly(P)
    w = deepcopy(x)
    g = Array{P}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem(pow_mod(w, prime, prime), f)(prime)
        g[k] = gcd(w - x, f, prime) 
        f = (f ÷ g[k])(prime)
    end


    #edge case for final factor
    f != one(P) && push!(g,f)
    
    return g
end

"""
Distinct degree split.
Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function dd_split(f :: PolynomialModP, d :: Int) :: Vector{PolynomialModP}
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand(PolynomialModP, degree = d, monic = true)
    n_power = ((f.prime)^d - 1) ÷ 2
    g = gcd(w^n_power - one(f), f)
    ḡ = f ÷ g
    return vcat(dd_split(g, d), dd_split(ḡ, d))
end

function dd_split(f :: P, d :: Int, prime :: Int) where P <: Polynomial
    f = mod(f,prime)
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand(P, degree = d, monic = true)
    w = mod(w,prime)
    n_power = (prime^d-1) ÷ 2
    g = gcd(pow_mod(w,n_power,prime) - one(P), f, prime)
    ḡ = (f ÷ g)(prime) # g\bar + [TAB]
    return vcat(dd_split(g, d, prime), dd_split(ḡ, d, prime) )
end
