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
include("polynomial_sparse128_test.jl")
include("polynomial_mod_p_test.jl")
include("factorization_test.jl")

####
# Execute unit tests for integers
###
# test_euclid_ints()
# test_ext_euclid_ints()

println("\n\n")
println("######################################")
println("# Polynomial mod p tests #")
println("######################################")
@time prod_test_poly_mod_p()
@time prod_derivative_test_poly_mod_p()
@time ext_euclid_test_poly_mod_p()
@time division_test_poly_mod_p()
@time factor_test_poly_mod_p()

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
@time pow_mod_poly()

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
@time pow_mod_poly_sparse()

####
# Execute unit tests for 128-bit sparse polynomials and factorization
####
println("\n\n")
println("######################################")
println("# 128-bit Sparse polynomial tests #")
println("######################################")
@time prod_test_poly_sparse_128()
@time prod_derivative_test_poly_sparse_128()
@time ext_euclid_test_poly_sparse_128()
@time division_test_poly_sparse_128()
@time factor_test_poly_sparse_128()
@time pow_mod_poly_sparse_128()
@time overflow_poly_sparse_128()

####
# Execute unit tests for polynomials mod p
####
println("\n\n")
println("######################################")
println("# Polynomial mod p tests #")
println("######################################")
@time prod_test_poly_mod_p()
@time prod_derivative_test_poly_mod_p()
@time ext_euclid_test_poly_mod_p()
@time division_test_poly_mod_p()
@time factor_test_poly_mod_p()
