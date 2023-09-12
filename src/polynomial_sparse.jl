include("dict_linked_list.jl")

"""
Sparse Polynomial type - store only the monomials with non-zero coefficient
"""
struct PolynomialSparse

    # DictLinkedList mapping degrees of terms to the term. A polynomial is the
    # zero polynomial if and only if DictLinkedList is empty
    terms :: DictLinkedList{Int, Term}
    
    #Inner constructor of 0 polynomial
    PolynomialSparse() = new(DictLinkedList{Int, Term}(isless))

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse(vt::Vector{Term})
        terms = DictLinkedList{Int, Term}(isless)
        for t in vt
            if iszero(t)
                continue
            elseif contains(terms, t.degree)
                replace!(terms, t.degree, t)
            else
                insert!(terms, t.degree, t)
            end
        end
        new(terms)
    end

    PolynomialSparse(dll :: DictLinkedList{Int, Term}) = new(dll)
end

"""
Construct a polynomial with a single term.
"""
PolynomialSparse(t::Term) = PolynomialSparse([t])

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial_sparse(p::Int) = PolynomialSparse([Term(1,p), Term(-1,0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial_sparse(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly(PolynomialSparse) = PolynomialSparse(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
one(p::PolynomialSparse) = one(PolynomialSparse)

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialSparse} ;
                degree::Int = -1,
                terms::Int = -1,
                max_coeff::Int = 100,
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)

    while true
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialSparse([Term(coeffs[i],degrees[i]) for i in 1:length(degrees)])
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialSparse) 
    if iszero(p)
        print(io,"0")
    else
        is_first = true
        n = length(p.terms)
        term_list = lowest_to_highest ? p.terms : p.terms |> collect |> reverse
        for t in term_list
            if is_first
                print(io, t.coeff < 0 ? "-" : "", t)
            else
                print(io, t.coeff < 0 ? " - " : " + ", t)
            end
            is_first = false
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial.
This implements the iteration interface.
"""
iterate(p::PolynomialSparse, state = p.terms.list.node) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse) = length(p.terms)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : last(p.terms)

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparse)::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::PolynomialSparse)::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
TODO
"""
evaluate(f::PolynomialSparse, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
push!(p::PolynomialSparse, t::Term) = insert!(p.terms, t.degree, t)

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialSparse)::Term 
    if iszero(p)
        return zero(Term)
    end
    term = leading(p)
    remove!(p.terms, term.degree)
    return term
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialSparse)::Bool = p.terms |> empty

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparse)::PolynomialSparse 
    der_p = PolynomialSparse()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return der_p
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialSparse) = p ÷ content(p)


"""
A square free polynomial.
"""
square_free(p::PolynomialSparse, prime::Int)::PolynomialSparse = (
    p ÷ gcd(p,derivative(p),prime))(prime)

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialSparse, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse = p1 + (-p2)


"""
Multiplication of polynomial and term.
"""
function *(t::Term, p1::PolynomialSparse)::PolynomialSparse
    if iszero(t) || iszero(p1)
        return PolynomialSparse()
    end
    p = PolynomialSparse()
    for term in p1
        push!(p, term * t)
    end
    return p
end
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod(f :: PolynomialSparse, p :: Int) :: PolynomialSparse
    f_out = PolynomialSparse(map(t -> mod(t, p), f.terms))
    filter!(f_out.terms, t -> !iszero(t))
    return f_out
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p :: PolynomialSparse, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end
