#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials modulo prime.
"""
function extended_euclid_alg(a :: PolynomialModP, b :: PolynomialModP)
    @assert a.prime == b.prime
    old_r, r = a, b
    old_s, s = one(a), zero(a)
    old_t, t = zero(a), one(a)

    while !iszero(r)
        q = first(divide(old_r, r))
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    end
    g, s, t = old_r, old_s, old_t

    @assert s*a + t*b - g == 0
    return g, s, t
end

function extended_euclid_alg(a :: P, b :: P, prime :: Int) where P <: Polynomial
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(P), zero(P)
    old_t, t = zero(P), one(P)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t
end


"""
The GCD of two polynomials modulo prime.
"""
function gcd(a :: PolynomialModP, b :: PolynomialModP)
    return extended_euclid_alg(a, b) |> first
end

function gcd(a::P, b::P, prime::Int) where P <: Polynomial
    return extended_euclid_alg(a, b, prime) |> first
end
