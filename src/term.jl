#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.
"""
struct Term{I}
    coeff :: I
    degree :: Int
    function Term{I}(coeff::I, degree::Int) where I <: Integer
        degree < 0 && error("Degree must be non-negative")
        iszero(coeff) ? new(coeff,0) : new(coeff,degree)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term{ResidueInt}}, prime :: Int) = Term{ResidueInt}(
    zero(ResidueInt, prime), 0)
zero(::Type{Term{I}}) where I <: Integer = Term{I}(I(0),0)

"""
Creates the unit term.
"""
one(::Type{Term{ResidueInt}}, prime :: Int) = Term{ResidueInt}(
    one(ResidueInt, prime),0)
one(::Type{Term{I}}) where I <: Integer = Term{I}(I(1),0)

###########
# Display #
###########

"""
Show a term.
"""
function show(io::IO, t::Term{I}) where I <: Integer
    if t.degree == 0
        print(io, abs(t.coeff))
    elseif abs(t.coeff) == 1
        print(io, "x", t.degree == 1 ? "" : int_to_superscript(t.degree))
    else
        print(io, abs(t.coeff), "x",
            t.degree == 1 ? "" : int_to_superscript(t.degree))
    end
end

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t :: Term{I}) where I <: Integer = iszero(t.coeff)

"""
Compare two terms. We pinky promise to only call this on Term{ResidueInt} when
we only care about ordering by degree
"""
function isless(t1::Term{ResidueInt},t2::Term{ResidueInt})
    return t1.degree < t2.degree
end

function isless(t1::Term{I},t2::Term{I}) where I <: Integer
    return t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  
end

"""
Evaluate a term at a point x.
"""
evaluate(t::Term{I}, x::T) where {T <: Number, I <: Integer} = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1 :: Term{I} ,t2 :: Term{I}) where I <: Integer
    @assert t1.degree == t2.degree
    Term{I}(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term{I},) where I <: Integer = Term{I}(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term{I}, t2::Term{I}) where I <: Integer = t1 + (-t2) 

"""
Multiply two terms.
"""
function *(t1::Term{I}, t2::Term{I}) where I <: Integer
    return Term{I}(t1.coeff * t2.coeff, t1.degree + t2.degree)
end

"""
Multiply a term by an integer. Could return 0
"""
function *(n :: J, t :: Term{I}) where {I <: Integer, J <: Integer}
    coefficient = n * t.coeff
    return iszero(coefficient) ? Term{I}(coefficient, 0) : Term{I}(coefficient, t.degree)
end


"""
Compute the symmetric mod of a term with an integer.
"""
function mod(t :: Term{ResidueInt}, p :: I) where I <: Integer
    return error("Cannot take mod in a field constructed using mod")
end
function mod(t::Term{I}, p::J) where {I <: Integer, J <: Integer}
    return Term{I}(mod(t.coeff, I(p)), t.degree)
end

"""
Compute the derivative of a term.
"""
derivative(t::Term{I}) where I <: Integer = Term{I}(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer if the coefficient ring is
not a field.
"""
function ÷(t1::Term{ResidueInt},t2::Term{ResidueInt})
    @assert t1.degree >= t2.degree
    return Term{ResidueInt}(t1.coeff ÷ t2.coeff, t1.degree - t2.degree)
end
function ÷(t1::Term{I},t2::Term{I}) where I <: Integer
    @assert t1.degree >= t2.degree
    f(p::J) where J <: Integer = Term{I}(mod((
        t1.coeff * int_inverse_mod(t2.coeff, I(p))), I(p)), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term{ResidueInt}, n::I) where I <: Integer = Term{ResidueInt}(t.coeff ÷ n, t.degree)
÷(t::Term{I}, n::J) where {I <: Integer, J <: Integer} = t ÷ Term{I}(I(n),0)
