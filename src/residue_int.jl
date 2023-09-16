"""
Representation of an integer modulo a prime
"""
struct ResidueInt <: Integer
    # Prime used to construct the field of residue classes
    prime :: Int
    # Value mod prime
    value :: Int

    function ResidueInt(value :: I, p :: J) where {I <: Integer, J <: Integer}
        return new(Int(p), Int(value % p))
    end
end

"""
Give the smallest strictly positive representative of the equivalence class
"""
positive_value(x :: ResidueInt) = x.value < 0 ? ResidueInt(x.value + x.prime, x.prime) : x

"""
Show the canonical representative of the int
"""
show(io :: IO, n :: ResidueInt) = print(io, "$(n.value)")


"""
Integer conversion methods
"""
Int(n :: ResidueInt) = n.value
Int128(n :: ResidueInt) = n.value
ResidueInt(n :: ResidueInt) = ResidueInt(n.value, n.prime)


"""
Unit and zero constructors
"""
zero(::Type{ResidueInt}, prime :: Int) = ResidueInt(0, prime)
one(::Type{ResidueInt}, prime :: Int) = ResidueInt(1, prime)
iszero(x :: ResidueInt) = x.value == 0

"""
ResidueInt acts kind of like a functor (loosely; it's not quite polymorphic
enough) so we may define fmap
"""
fmap(f :: Function, x :: ResidueInt) = ResidueInt(f(x.value), x.prime)

"""
Basic arithmetic operations
"""
function +(x :: ResidueInt, y :: ResidueInt)
    @assert x.prime == y.prime
    return ResidueInt(x.value + y.value, x.prime)
end

function -(x :: ResidueInt, y :: ResidueInt)
    @assert x.prime == y.prime
    return ResidueInt(x.value - y.value, x.prime)
end

function *(x :: ResidueInt, y :: ResidueInt)
    @assert x.prime == y.prime
    return ResidueInt(x.value * y.value, x.prime)
end

function inverse(x :: ResidueInt)
    @assert !iszero(x)
    g, t, _ = ext_euclid_alg(x.value, x.prime)
    @assert g == 1
    return ResidueInt(t, x.prime)
end

function ÷(x :: ResidueInt, y :: ResidueInt)
    @assert x.prime == y.prime
    return ResidueInt(x.value * inverse(y), y.prime)
end

-(x :: ResidueInt) = fmap(-, x)

-(x :: ResidueInt, y :: ResidueInt) = x + (-y)


"""
Only to be use for printing
"""
abs(x :: ResidueInt) = x.value < 0 ? -x : x

"""
Equality of integers up to modulo a prime
"""
==(x :: ResidueInt, y :: ResidueInt) = iszero(x - y)
==(x :: ResidueInt, y :: I) where I <: Integer = x == ResidueInt(y, x.prime)
==(y :: I, x :: ResidueInt) where I <: Integer = x == y

"""
Arithmetic operations with integers
"""
+(x :: ResidueInt, y :: I) where I <: Integer = x + ResidueInt(y, x.prime)
+(y :: I, x :: ResidueInt) where I <: Integer = x + ResidueInt(y, x.prime)

*(x :: ResidueInt, y :: I) where I <: Integer = x * ResidueInt(y, x.prime)
*(y :: I, x :: ResidueInt) where I <: Integer = x * ResidueInt(y, x.prime)

-(x :: ResidueInt, y :: I) where I <: Integer = x - ResidueInt(y, x.prime)
-(y :: I, x :: ResidueInt) where I <: Integer = ResidueInt(y, x.prime) - x

÷(x :: ResidueInt, y :: I) where I <: Integer = x ÷ ResidueInt(y, x.prime)

function ^(x :: ResidueInt, y :: I) where I <: Integer
    iszero(x) && return deepcopy(x)
    y < 0 && return inverse(x) ^ (-y)
    ret = one(ResidueInt, x.prime)
    d = y
    while d > 0
        ret *= x
        d -= 1
        ret == 1 && return x^(y % (y - d)) # x^(y - d) = 1
    end
    return ret
end
