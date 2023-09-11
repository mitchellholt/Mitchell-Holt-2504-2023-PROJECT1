#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::Polynomial, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end

function +(p::PolynomialSparse, t::Term)
    p_out = deepcopy(p)
    if iszero(t)
        return p_out
    end
    if degree(p_out) < t.degree || !contains(p_out.terms, t.degree)
        push!(p_out, t)
        return p_out
    else
        replace!(t.degree, lookup(p_out.terms, t.degree) + t)
        # Make sure we haven't just added a zero term
        if iszero(lookup(p_out.terms, t.degree))
            remove!(p_out.terms, t.degree)
        end
        return p_out
    end
end

+(t::Term, p::Union{Polynomial, PolynomialSparse}) = p + t

"""
Add two polynomials.
"""
function +(p1 :: P, p2 :: P) :: P where P <: Union{Polynomial, PolynomialSparse}
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Union{Polynomial, PolynomialSparse}, n::Int) = p + Term(n,0)
+(n::Int, p::Union{Polynomial, PolynomialSparse}) = p + Term(n,0)
