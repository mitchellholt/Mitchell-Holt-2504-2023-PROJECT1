#############################################################################
#############################################################################
#
# This file implement several basic algorithms for integers. 
#                                                                               
#############################################################################
#############################################################################

global use_unicode = true

"""
The quotient between two numbers
"""
function quo(a::T,b::T) where T <: Integer
    a < 0 && return -quo(-a,b)
    a < b && return 0
    return 1 + quo(a-b, b)
end

"""
Symmetric mod over an odd number
"""
function smod(n :: I, p :: I) where I <: Integer
    @assert p % 2 == 1
    k = div(p, 2)
    r = mod(n, p)
    return r <= k ? r : r - p
end

"""
Euclid's algorithm.
"""
function euclid_alg(a, b; rem_function = %)
    b == 0 && return a
    return euclid_alg(b, rem_function(a,b); rem_function = rem_function)
end

"""
Euclid's algorithm on multiple arguments.
"""
euclid_alg(a...) = foldl((a,b)->euclid_alg(a,b), a; init = 0)

"""
Euclid's algorithm on a vector.
"""
euclid_alg(a::Vector{T}) where T <: Integer = euclid_alg(a...)


"""
The extended Euclidean algorithm.
"""
function ext_euclid_alg(a, b, rem_function = %, div_function = ÷)
    a == 0 && return b, 0, 1
    g, t, s = ext_euclid_alg(rem_function(b,a), a, rem_function, div_function)
    s = s - div_function(b,a)*t
    @assert g == a*s + b*t
    return g, s, t
end

"""
Raise an integer to a power modulo a prime
"""
function pow_mod(k :: I, n :: Int, prime :: J) where {I <: Integer, J <: Integer}
    res = 1
    while n > 0
        res = (res * k) % I(prime)
        n -= 1
    end
    return res
end

"""
Display the result of the extended Euclidean algorithm.
"""
pretty_print_egcd((a,b),(g,s,t)) = println("$a × $s + $b × $t = $g") #\times + [TAB]

"""
Integer inverse symmetric mod
"""
function int_inverse_mod(a :: I, m :: I) where I <: Integer
    if mod(a, m) == 0
        error("Can't find inverse of $a mod $m because $m divides $a") 
    end
    return mod(ext_euclid_alg(a,m)[2],m)
end

"""
Convert a single base 10 digit to a superscript unicode character
"""
function digit_to_superscript(n :: Int) :: Char
    # Literal for zero
    zero_code = 0x2070
    
    # cringe edge cases
    if n == 1
        return Char(0x00B9)
    elseif n == 2
        return Char(0x00B2)
    elseif n == 3
        return Char(0x00B3)
    else
        # based normal case
        return Char(zero_code + n)
    end
end

"""
Convert a positive integer to a string of superscipt unicode characters.
"""
function int_to_superscript(k :: Int) :: String
    if use_unicode
        return k == 0 ? "" : int_to_superscript(div(k, 10)) * digit_to_superscript(mod(k, 10))
    else
        return k == 0 ? "^" : int_to_superscript(div(k, 10)) * Char(Int('0') + mod(k, 10))
    end
end
