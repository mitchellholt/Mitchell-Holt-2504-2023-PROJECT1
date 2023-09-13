#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

include("src/general_alg.jl")
include("src/term.jl")

# Abstract polynomial type
abstract type Polynomial end

include("src/polynomial_dense.jl")
include("src/polynomial_sparse.jl")

# Type aliases
const PolynomialSparse = PolynomialSparse_{Int}
const PolynomialSparse128 = PolynomialSparse_{Int128}

include("src/basic_polynomial_operations/polynomial_addition.jl")
include("src/basic_polynomial_operations/polynomial_multiplication.jl")
include("src/basic_polynomial_operations/polynomial_division.jl")
include("src/basic_polynomial_operations/polynomial_gcd.jl")
include("src/polynomial_factorization/factor.jl")

nothing
