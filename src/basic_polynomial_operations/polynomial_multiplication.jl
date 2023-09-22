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
function crt(a :: PolynomialModP, b :: PolynomialModP)
    @assert euclid_alg(a.prime, b.prime) == 1
    x = x_poly(PolynomialModP, a.prime * b.prime)
    c = zero(x)
    for k in max(degree(a), degree(b)):(-1):0
        # Coefficient of x^k in a
        ak = contains(a.terms, k) ? lookup(a.terms, k).coeff : zero(ResidueInt, a.prime)
        # Coefficient of x^k in b
        bk = contains(b.terms, k) ? lookup(b.terms, k).coeff : zero(ResidueInt, b.prime)
        # Should be ResidueInt with "prime" a.prime * b.prime
        ck = crt(ak, bk)
        c += ck * x^k
    end
    return c
end


"""
Use polynomial CRT to multiply two polynomials
"""
function multiply(a :: PolynomialSparse128, b :: PolynomialSparse128)
    # Create an upper bound for the prime to use
    height_a = maximum(abs(t.coeff) for t in a)
    height_b = maximum(abs(t.coeff) for t in b)
    upper = 2 * height_a * height_b * min(degree(a) + 1, degree(b) + 1)

    p = Int128(3)
    prime_product = p
    c = PolynomialModP(a, p) * PolynomialModP(b, p) # not mod prime, but its ok for here
    while prime_product < upper
        p = next_prime(p)
        c = crt(c, PolynomialModP(a, p) * PolynomialModP(b, p))
        prime_product *= p
    end
    return PolynomialSparse128(c)
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
Power of a polynomial using repeated squaring
"""
function ^(p :: P, n :: I) where {P <: Polynomial, I <: Integer}
    n < 0 && error("No negative power")
    m = UInt128(n)
    ans, w = one(p), p
    while n > 0
        if n & 0x01 != 0x00
            ans *= w
        end
        w *= w
        n >>= 1
    end
    return ans
end
