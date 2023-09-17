#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num :: PolynomialDense, den :: PolynomialDense)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialDense()
        prev_degree = degree(f)
        while degree(f) >= degree(g) 
            h = PolynomialDense((leading(f) ÷ leading(g))(p))  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

function divide(num :: PolynomialModP, den :: PolynomialModP)
    f = deepcopy(num)
    g = den
    iszero(g) && throw(DivideError())
    q = PolynomialModP()
    prev_degree = degree(f)
    while degree(f) >= degree(g) 
        h = PolynomialModP(leading(f) ÷ leading(g))
        f = f - h*g
        q += h
        prev_degree == degree(f) && break
        prev_degree = degree(f)
    end
    @assert iszero(num  - (q*g + f))
    return q, f
end

function divide(num :: PolynomialSparse_{I}, den :: PolynomialSparse_{I}) where I <: Integer
    function division_function(prime::J) where J <: Integer
        p = I(prime)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparse_{I}()
        prev_degree = degree(f)
        while degree(f) >= degree(g) 
            h = PolynomialSparse_{I}((leading(f) ÷ leading(g))(p))  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end
"""
The quotient from polynomial division. Returns a function of an integer.
"""
function ÷(num :: PolynomialModP, den :: PolynomialModP)
    first(divide(num, den))
end

function ÷(num :: P, den :: P) where P <: Polynomial
    return p -> first(divide(num,den)(p))
end

"""
The remainder from polynomial division. Returns a function of an integer.
"""
function rem(num :: PolynomialModP, den :: PolynomialModP)
    last(divide(num,den))
 end

function rem(num :: P, den :: P) where P <: Polynomial
    return p -> last(divide(num,den)(p))
 end
