"""
Sparse Polynomial type with integer coeffiecients modulo a prime
"""
struct PolynomialModP <: Polynomial

    # DictLinkedList mapping degrees of terms to the term. A polynomial is the
    # zero polynomial if and only if DictLinkedList is empty
    terms :: DictLinkedList{Int, Term{ResidueInt}}
    prime :: Int
    
    #Inner constructor of 0 polynomial
    PolynomialModP(prime :: Int) = new(
        DictLinkedList{Int, Term{ResidueInt}}(isless), prime)

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialModP(vt::Vector{Term{ResidueInt}}, prime :: Int)
        terms = DictLinkedList{Int, Term{ResidueInt}}(isless)
        for t in vt
            if iszero(t)
                continue
            elseif contains(terms, t.degree)
                replace!(terms, t.degree, t)
            else
                insert!(terms, t.degree, t)
            end
        end
        new(terms, prime)
    end

    function PolynomialModP(vt::Vector{Term{I}}, prime :: Int) where I <: Integer
        terms = DictLinkedList{Int, Term{ResidueInt}}(isless)
        for t in map(term -> Term{ResidueInt}(
                ResidueInt(term.coeff, prime), term.degree), vt)
            if iszero(t)
                continue
            elseif contains(terms, t.degree)
                replace!(terms, t.degree, t)
            else
                insert!(terms, t.degree, t)
            end
        end
        new(terms, prime)
    end

    PolynomialModP(dll :: DictLinkedList{Int, Term{ResidueInt}}, prime :: Int) = new(dll, prime)
end

"""
Construct a polynomial with a single term.
"""
PolynomialModP(t::Term{ResidueInt}) = PolynomialModP([t], t.coeff.prime)

"""
Construct a polynomial of the form x^p-x.
"""
function cyclotonic_polynomial(::Type{PolynomialModP}, prime :: Int)
    return PolynomialModP([
            Term{ResidueInt}(ResidueInt(1, prime), prime),
            Term{ResidueInt}(ResidueInt(-1, prime), 0)],
        prime)
end

"""
Construct a polynomial of the form x-n.
"""
function linear_monic_polynomial(::Type{PolynomialModP}, n::I, prime :: Int) where I <: Integer
    return PolynomialModP([
            Term{ResidueInt}(ResidueInt(1, prime), 1),
            Term{ResidueInt}(ResidueInt(-n, prime), 0)],
        prime)
end

"""
Construct a polynomial of the form x.
"""
function x_poly(::Type{PolynomialModP}, p :: Int)
    return PolynomialModP(Term{ResidueInt}(ResidueInt(1, p), 1))
end

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialModP}, prime) = PolynomialModP(prime)
zero(p :: PolynomialModP) = PolynomialModP(p.prime)

"""
Creates the unit polynomial.
"""
function one(::Type{PolynomialModP}, prime :: Integer)
    PolynomialModP(one(Term{ResidueInt}, prime))
end
one(p::PolynomialModP) = one(PolynomialModP, p.prime)

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialModP}, prime :: Int;
              degree::Int = -1,
              terms::Int = -1,
              mean_degree::Float64 = 5.0,
              prob_term::Float64  = 0.7,
              monic = false,
              condition = (p)->true)

    while true
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = [ResidueInt(n, prime) for n in rand(1:(prime - 1),_terms+1)]
        monic && (coeffs[end] = one(ResidueInt, prime))
        p = PolynomialModP(
            [Term{ResidueInt}(coeffs[i],degrees[i]) for i in 1:length(degrees)],
            prime)
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialModP)
    if iszero(p)
        print(io,"0")
    else
        is_first = true
        n = length(p.terms)
        term_list = lowest_to_highest ? p.terms : p.terms |> collect |> reverse
        for t in term_list
            if is_first
                print(io, t.coeff.value < 0 ? "-" : "", t)
            else
                print(io, t.coeff.value < 0 ? " - " : " + ", t)
            end
            is_first = false
        end
        print(io, "  (mod $(p.prime))")
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial.
This implements the iteration interface.
"""
function iterate(p::PolynomialModP, state = p.terms.list.node)
    return iterate(p.terms, state)
end

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialModP) = length(p.terms)

"""
The leading term of the polynomial.
"""
function leading(p::PolynomialModP)
    return isempty(p.terms) ? zero(Term{ResidueInt}, p.prime) : last(p.terms)
end

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialModP) = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::PolynomialModP) = iszero(p) ? 0 : leading(p).degree 

"""
Evaluate the polynomial at a point `x`. Can only do this over integer-like types
"""
function evaluate(f :: PolynomialModP, x :: I) where I <: Integer
    return sum(evaluate(t,x) for t in f)
end

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
push!(p::PolynomialModP, t::Term{ResidueInt}) = insert!(p.terms, t.degree, t)

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialModP)
    if iszero(p)
        return zero(Term{ResidueInt}, p.prime)
    end
    term = leading(p)
    remove!(p.terms, term.degree)
    return term
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialModP) = p.terms |> empty

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
function -(p::PolynomialModP)
    PolynomialModP(map(-, p.terms), p.prime)
end

"""
Create a new polynomial which is the derivative of the polynomial. This doesn't
really make sense for a polynomial over a field of residue classes of the
integers modulo a prime, but doing a mod operation before or after taking a
derivative makes no difference.
"""
function derivative(p::PolynomialModP)
    der_p = PolynomialModP()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return der_p
end


"""
A square free polynomial.
"""
square_free(p::PolynomialModP) = p ÷ gcd(p, derivative(p))

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialModP, p2::PolynomialModP) = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialModP, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
function -(p1::PolynomialModP, p2::PolynomialModP)
    return p1 + (-p2)
end


"""
Multiplication of polynomial and term.
"""
function *(t::Term{ResidueInt}, p1::PolynomialModP)
    p = PolynomialModP(p1.prime)
    (iszero(t) || iszero(p1)) && return p
    for term in p1
        push!(p, term * t)
    end
    return p
end
*(p1::PolynomialModP, t::Term{ResidueInt}) = t*p1

"""
Multiplication of polynomial and an integer.
"""
function *(n::I, p1::PolynomialModP) where I <: Integer
    p = PolynomialModP(p1.prime)
    i = ResidueInt(n, p1.prime)
    (iszero(i) || iszero(p1)) && return p
    for term in p1
        t = i * term
        push!(p, t)
    end
    return p
end
*(p::PolynomialModP, n::I) where I <: Integer = n*p

"""
Integer division of a polynomial by an integer.
Warning this may not make sense if n does not divide all the coefficients of p.
Return a function if we are not working over a field
"""
function ÷(p::PolynomialModP, n::I) where I <: Integer
    return map(pt -> pt ÷ n, p.terms)
end
