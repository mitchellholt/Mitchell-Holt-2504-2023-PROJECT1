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
function divide(num :: P, den :: P) where P <: Union{PolynomialDense, PolynomialSparse}
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = P()
        prev_degree = degree(f)
        while degree(f) >= degree(g) 
            h = P((leading(f) ÷ leading(g))(p))  #syzergy 
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
function ÷(num :: P, den :: P) where P <: Union{PolynomialDense, PolynomialSparse}
    return (p :: Int) -> divide(num,den)(p) |> first
end

"""
The remainder from polynomial division. Returns a function of an integer.
"""
function rem(num :: P, den :: P) where P <: Union{PolynomialDense, PolynomialSparse}
     return (p :: Int) -> divide(num,den)(p) |> last
 end
