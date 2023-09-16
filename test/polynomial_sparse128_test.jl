#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations with the sparse
# polynomial data type
#                                                                               
#############################################################################
#############################################################################

"""
Test product of polynomials.
"""
function prod_test_poly_sparse_128(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end
    for _ in 1:N
        p_base = one(PolynomialSparse128)
        for _ in 1:N_prods
            p = rand(PolynomialSparse128)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_sparse_128 - PASSED")
end

"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly_sparse_128(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly_sparse_128 - PASSED")
end


"""
Test division of polynomials modulo p.
"""
function division_test_poly_sparse_128(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p_prod = p1*p2
        q, r = PolynomialSparse128(), PolynomialSparse128()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_poly_sparse_128 - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly_sparse_128(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly_sparse_128 - PASSED")
end


"""
Super sparse polynomials

"""
function pow_mod_poly_sparse_128(;prime::Int = 19, N::Int = 30, seed::Int = 0)
    Random.seed!(seed)
    print("pow_mod_poly_sparse_128\t")
    for k in 1:N
        f = rand(PolynomialSparse128; prob_term = 0.01, mean_degree = Float64(120))
        @assert leading(pow_mod(f, k, prime)).coeff == pow_mod(leading(f).coeff, k, prime)
        print(".")
    end
    println("\tPASSED")
end

"""
Test using integer arithmetic which would overflow in the PolynomialSparse type
"""
function overflow_poly_sparse_128()
    x = x_poly(PolynomialSparse128)
    f = Int128(2^50) * x
    @assert leading(f^2).coeff == (leading(f).coeff)^2
end
