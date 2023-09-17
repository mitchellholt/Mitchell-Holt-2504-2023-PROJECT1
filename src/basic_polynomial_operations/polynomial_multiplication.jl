#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

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
