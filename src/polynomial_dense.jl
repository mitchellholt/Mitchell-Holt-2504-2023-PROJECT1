#############################################################################
#############################################################################
#
# This file defines the polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# PolynomialDense type and construction #
####################################

"""
Based global variable for deciding the order in which polynomial terms are
printed
"""
global lowest_to_highest :: Bool = false

"""
A PolynomialDense type - designed to be for polynomials with integer coefficients.
"""
struct PolynomialDense <: Polynomial

    #A zero packed vector of terms
    #Term{Int}s are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    #Note: at positions where the coefficient is 0, the power of the term is also 0 (this is how the Term{Int} type is designed)
    terms::Vector{Term{Int}}   
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term{Int})])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term{Int}})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term{Int})]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term{Int}) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new(terms)
    end
end

"""
This function maintains the invariant of the PolynomialDense type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::PolynomialDense)::PolynomialDense
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

"""
Construct a polynomial with a single term.
"""
PolynomialDense(t::Term{Int}) = PolynomialDense([t])

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int) = PolynomialDense([Term{Int}(1,p), Term{Int}(-1,0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int) = PolynomialDense([Term{Int}(1,1), Term{Int}(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly(PolynomialDense) = PolynomialDense(Term{Int}(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialDense})::PolynomialDense = PolynomialDense()
zero(p :: PolynomialDense)::PolynomialDense = PolynomialDense()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialDense})::PolynomialDense = PolynomialDense(one(Term{Int}))
one(p::PolynomialDense) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialDense} ; 
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
        p = PolynomialDense( [Term{Int}(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialDense) 
    if iszero(p)
        print(io,"0")
    else
        is_first = true
        n = length(p.terms)
        term_list = lowest_to_highest ? p.terms : reverse(p.terms)
        for t in term_list
            if iszero(t)
                continue
            elseif is_first
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
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialDense, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialDense) = length(p.terms) 

"""
The leading term of the polynomial.
"""
leading(p::PolynomialDense)::Term{Int} = isempty(p.terms) ? zero(Term{Int}) : last(p.terms) 

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::PolynomialDense)::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::PolynomialDense)::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialDense)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialDense, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialDense, t::Term{Int}) 
    if t.degree <= degree(p)
        p.terms[t.degree + 1] = t
    else
        append!(p.terms, zeros(Term{Int}, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialDense)::Term{Int} 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term{Int}))
    end

    return popped_term
end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialDense)::Bool = p.terms == [Term{Int}(0,0)]

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialDense) = PolynomialDense(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialDense)::PolynomialDense 
    der_p = PolynomialDense()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return trim!(der_p)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialDense) = p ÷ content(p)


"""
A square free polynomial.
"""
square_free(p::PolynomialDense, prime::Int)::PolynomialDense = (p ÷ gcd(p,derivative(p),prime))(prime)

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialDense, p2::PolynomialDense)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialDense, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense = p1 + (-p2)


"""
Multiplication of polynomial and term.
"""
*(t::Term{Int}, p1::PolynomialDense)::PolynomialDense = iszero(t) ? PolynomialDense() : PolynomialDense(map((pt)->t*pt, p1.terms))
*(p1::PolynomialDense, t::Term{Int})::PolynomialDense = t*p1

"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialDense)::PolynomialDense = p*Term{Int}(n,0)
*(p::PolynomialDense, n::Int)::PolynomialDense = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialDense, n::Int) = (prime)->PolynomialDense(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::PolynomialDense, p::Int)::PolynomialDense
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return trim!(f_out)
        
    # p_out = PolynomialDense()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::PolynomialDense, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end
