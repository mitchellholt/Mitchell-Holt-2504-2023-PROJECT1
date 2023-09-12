#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")
include("integers_test.jl")
include("polynomials_test.jl")
include("polynomial_sparse_test.jl")
include("factorization_test.jl")

println("Starting tests")

####
# Execute unit tests for integers
###
# test_euclid_ints()
# test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
@time prod_test_poly()
prod_derivative_test_poly()
ext_euclid_test_poly()
division_test_poly()

####
# Execute unit tests for sparse polynomials
####
@time prod_test_poly_sparse()
prod_derivative_test_poly_sparse()
ext_euclid_test_poly_sparse()
division_test_poly_sparse()

####
# Execute unit tests for polynomial factorization
####
@time factor_test_poly()
@time factor_test_poly_sparse()
