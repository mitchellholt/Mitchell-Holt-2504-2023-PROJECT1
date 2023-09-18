#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Chinese remainder theorem on two (mod p) polynomials. Construct a polynomial c
such that c = a (mod p) and c = b (mod q), where p = a.prime, q = b.prime.
Confusingly enough, a.prime and b.prime need not actually be prime, we only
require that they have gcd 1
"""
function crt(p1 :: PolynomialModP, p2 :: PolynomialModP)
    @assert euclid_alg(p1.prime, p2.prime) == 1
    x = x_poly(PolynomialSparse128)
    c = zero(PolynomialSparse128)
    for k in max(degree(a), degree(b)):(-1):0
        # Coefficient of x^k in a
        ak = contains(a.terms, k) ? a.terms[k].data.coeff : zero(ResidueInt, p1.prime)
        # Coefficient of x^k in b
        bk = contains(b.terms, k) ? b.terms[k].data.coeff : zero(ResidueInt, p2.prime)
        # Should be Int128
        ck = crt(ak, bk)
        c = c + ck * x
    end
end


"""
Use polynomial CRT to multiply two polynomials
"""
function multiply(a :: PolynomialSparse128, b :: PolynomialSparse128)
    # Create an upper bound for the prime to use
    height_a = maximum(t.coeff for t in a)
    height_b = maximum(t.coeff for t in b)
    n, m = degree(a), degree(b)
    upper = 2 * height_a * height_b * min(n, m)

    p = 3
    prime_product = Int128(p)
    c = PolynomialModP(a * b, prime_product) # not mod prime, but its ok for here
    while prime_product < upper
        p = next_prime(p)
        new_c = PolynomialModP(a * b, p)
        c = crt(c, new_c)
        prime_product *= p
    end

    return Polynomial128(c, prime_product)
end


"""
Multiply two polynomials.
"""
function *(p1 :: P, p2 :: P) :: P where P <: Polynomial
    p_out = zero(p1)
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial.
"""
function ^(p :: P, n :: Int) where P <: Polynomial
    n < 0 && error("No negative power")
    iszero(p) && return zero(p)
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end
