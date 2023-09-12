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
function +(p::PolynomialDense, t::Term)
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
        term = t + lookup(p.terms, t.degree)
        # Make sure we haven't just added a zero term
        if !iszero(term)
            replace!(p_out.terms, term.degree)
        end
        return p_out
    end
end

+(t::Term, p::Union{PolynomialDense, PolynomialSparse}) = p + t

"""
Add two polynomials.
"""
function +(p1 :: PolynomialDense, p2 :: PolynomialDense) :: PolynomialDense
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

# TODO speed this up using insert at end
function +(p1 :: PolynomialSparse, p2 :: PolynomialSparse) :: PolynomialSparse
    p1_ = deepcopy(p1) # expensive?
    p2_ = deepcopy(p2)
    p = PolynomialSparse()
    while !iszero(p1_) || !iszero(p2_)
        if (degree(p1_) == degree(p2_))
            t = pop!(p1_) + pop!(p2_)
            iszero(t) ? continue : push!(p, t)
        elseif degree(p1_) > degree(p2_)
            push!(p, pop!(p1_))
        else
            push!(p, pop!(p2_))
        end
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Union{PolynomialDense, PolynomialSparse}, n::Int) = p + Term(n,0)
+(n::Int, p::Union{PolynomialDense, PolynomialSparse}) = p + Term(n,0)
