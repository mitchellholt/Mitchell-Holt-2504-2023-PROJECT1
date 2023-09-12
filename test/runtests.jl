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

####
# Execute unit tests for integers
###
# test_euclid_ints()
# test_ext_euclid_ints()

####
# Execute unit tests for dense polynomials and factorization
####
println("######################################")
println("# Dense polynomial tests #")
println("######################################")
@time prod_test_poly()
@time prod_derivative_test_poly()
@time ext_euclid_test_poly()
@time division_test_poly()
@time factor_test_poly()

####
# Execute unit tests for sparse polynomials and factorization
####
println("\n\n")
println("######################################")
println("# Sparse polynomial tests #")
println("######################################")
@time prod_test_poly_sparse()
@time prod_derivative_test_poly_sparse()
@time ext_euclid_test_poly_sparse()
@time division_test_poly_sparse()
@time factor_test_poly_sparse()
