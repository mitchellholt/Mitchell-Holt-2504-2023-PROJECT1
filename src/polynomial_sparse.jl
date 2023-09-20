"""
Sparse Polynomial type - store only the monomials with non-zero coefficient
"""
struct PolynomialSparse_{I} <: Polynomial

    # DictLinkedList mapping degrees of terms to the term. A polynomial is the
    # zero polynomial if and only if DictLinkedList is empty
    terms :: DictLinkedList{Int, Term{I}}
    
    #Inner constructor of 0 polynomial
    PolynomialSparse_{I}() where I <: Integer = new(DictLinkedList{Int, Term{I}}(isless))

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse_{I}(vt::Vector{Term{I}}) where I <: Integer
        terms = DictLinkedList{Int, Term{I}}(isless)
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

    PolynomialSparse_{I}(dll :: DictLinkedList{Int, Term{I}}) where I <: Integer = new(dll)
end

"""
Construct a polynomial with a single term.
"""
PolynomialSparse_{I}(t::Term{I}) where I <: Integer = PolynomialSparse_{I}([t])

"""
Construct a polynomial of the form x^p-x.
"""
function cyclotonic_polynomial(::Type{PolynomialSparse_{I}}, p::Int) where I <: Integer
    return PolynomialSparse_{I}([Term{I}(1,p), Term{I}(-1,0)])
end

"""
Construct a polynomial of the form x-n.
"""
function linear_monic_polynomial(::Type{PolynomialSparse_{I}}, n::J) where {I <: Integer, J <: Integer}
    return PolynomialSparse_{I}([Term{I}(1,1), Term{I}(I(-n),0)])
end

"""
Construct a polynomial of the form x.
"""
x_poly(::Type{PolynomialSparse_{I}}) where I <: Integer = PolynomialSparse_{I}(Term{I}(I(1),1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialSparse_{I}}) where I <: Integer = PolynomialSparse_{I}()
zero(p :: PolynomialSparse_{I}) where I <: Integer = PolynomialSparse_{I}()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialSparse_{I}}) where I <: Integer = PolynomialSparse_{I}(one(Term{I}))
one(p::PolynomialSparse_{I}) where I <: Integer = one(PolynomialSparse_{I})

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialSparse_{I}};
              degree::Int = -1,
              terms::Int = -1,
              max_coeff::I = I(100),
              mean_degree::Float64 = 5.0,
              prob_term::Float64  = 0.7,
              monic = false,
              condition = (p)->true) where I <: Integer

    while true
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:I(max_coeff),_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialSparse_{I}([Term{I}(coeffs[i],degrees[i]) for i in 1:length(degrees)])
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialSparse_{I}) where I <: Integer
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
function iterate(p::PolynomialSparse_{I}, state = p.terms.list.node) where I <: Integer
    return iterate(p.terms, state)
end

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse_{I}) where I <: Integer = length(p.terms)

"""
The leading term of the polynomial.
"""
function leading(p::PolynomialSparse_{I}) where I <: Integer
    return isempty(p.terms) ? zero(Term{I}) : last(p.terms)
end

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparse_{I}) where I <: Integer = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::PolynomialSparse_{I}) where I <: Integer = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse_{I}) where I <: Integer = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
function evaluate(f::PolynomialSparse_{I}, x::T) where {I <: Integer, T <: Number}
    return sum(evaluate(t,x) for t in f)
end

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
push!(p::PolynomialSparse_{I}, t::Term{I}) where I <: Integer = insert!(p.terms, t.degree, t)

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialSparse_{I}) where I <: Integer
    if iszero(p)
        return zero(Term{I})
    end
    term = leading(p)
    remove!(p.terms, term.degree)
    return term
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialSparse_{I}) where I <: Integer = p.terms |> empty

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
function -(p::PolynomialSparse_{I}) where I <: Integer
    PolynomialSparse_{I}(map((pt)->-pt, p.terms))
end

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparse_{I}) where I <: Integer
    der_p = PolynomialSparse_{I}()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return der_p
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialSparse_{I}) where I <: Integer = p ÷ content(p)


"""
A square free polynomial.
"""
function square_free(p::PolynomialSparse_{I}, prime::J) where {I <: Integer, J <: Integer}
    (p ÷ gcd(p,derivative(p),I(prime)))(I(prime))
end

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
function ==(p1::PolynomialSparse_{I}, p2::PolynomialSparse_{I}) where I <: Integer
    p1.terms == p2.terms
end


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialSparse_{I}, n::T) where {I <: Integer, T <: Real} = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
function -(p1::PolynomialSparse_{I}, p2::PolynomialSparse_{I}) where I <: Integer
    return p1 + (-p2)
end


"""
Multiplication of polynomial and term.
"""
function *(t::Term{I}, p1::PolynomialSparse_{I}) where I <: Integer
    if iszero(t) || iszero(p1)
        return PolynomialSparse_{I}()
    end
    p = PolynomialSparse_{I}()
    for term in p1
        push!(p, term * t)
    end
    return p
end
*(p1::PolynomialSparse_{I}, t::Term{I}) where I <: Integer = t*p1

"""
Multiplication of polynomial and an integer.
"""
*(n::J, p::PolynomialSparse_{I}) where {I <: Integer, J <: Integer} = p*Term{I}(I(n),0)
*(p::PolynomialSparse_{I}, n::J) where {I <: Integer, J <: Integer} = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
function ÷(p::PolynomialSparse_{I}, n::J) where {I <: Integer, J <: Integer}
    return (prime)->PolynomialSparse_{I}(map((pt)->((pt ÷ I(n))(I(prime))), p.terms))
end

"""
Take the mod of a polynomial with an integer.
"""
function mod(f :: PolynomialSparse_{I}, p :: J) where {I <: Integer, J <: Integer}
    f_out = PolynomialSparse_{I}(map(t -> mod(t, I(p)), f.terms))
    filter!(f_out.terms, t -> !iszero(t))
    return f_out
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p :: PolynomialSparse_{I},
        n :: Int, prime :: J) where {I <: Integer, J <: Integer}
    n < 0 && error("No negative power")
    return PolynomialSparse_{I}(PolynomialModP(p, prime)^n, mod)
end
