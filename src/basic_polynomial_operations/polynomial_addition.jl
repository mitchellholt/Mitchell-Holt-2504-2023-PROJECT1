#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a dense polynomial and a term.
"""
function +(p::PolynomialDense, t::Term{Int})
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

"""
Add a sparse polynomial and a term.
"""
function +(p::PolynomialSparse_{I}, t::Term{I}) where I <: Integer
    p_out = deepcopy(p)
    iszero(t) && p_out
    if !contains(p_out.terms, t.degree)
        push!(p_out, t)
        return p_out
    else
        term = t + lookup(p.terms, t.degree)
        # Make sure we haven't just added a zero term
        if !iszero(term)
            replace!(p_out.terms, term.degree, term)
        end
        return p_out
    end
end

function +(p::PolynomialModP, t::Term{ResidueInt})
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
            replace!(p_out.terms, term.degree, term)
        end
        return p_out
    end
end

+(t::Term{I}, p::PolynomialSparse_{I}) where I <: Integer = p + t
+(t::Term{ResidueInt}, p::PolynomialModP) = p + t
+(t::Term{Int}, p::PolynomialDense) = p + t

"""
Add two dense polynomials.
"""
function +(p1 :: PolynomialDense, p2 :: PolynomialDense) :: PolynomialDense
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add two sparse polynomials. We exclusively use heap-like operations for efficiency
"""
function +(p1 :: PolynomialSparse_{I}, p2 :: PolynomialSparse_{I}) where I <: Integer
    p1_ = deepcopy(p1)
    p2_ = deepcopy(p2)
    p = PolynomialSparse_{I}()
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

function +(p1 :: PolynomialModP, p2 :: PolynomialModP) where I <: Integer
    @assert p1.prime == p2.prime
    p1_ = deepcopy(p1)
    p2_ = deepcopy(p2)
    p = PolynomialModP(p1.prime)
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
function +(p :: PolynomialModP, n :: I) where I <: Integer
    return p + (n * one(Term{ResidueInt}, p.prime))
end
+(p :: PolynomialSparse_{I}, n :: I) where I <: Integer = p + Term{I}(n,0)
+(n :: I, p :: PolynomialSparse_{I}) where I <: Integer = p + Term{I}(n,0)
+(p :: PolynomialDense, n :: Int) = p + Term{Int}(n, 0)
+(n :: Int, p :: PolynomialDense) = p + Term{Int}(n, 0)
