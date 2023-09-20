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
function prod_test_poly_mod_p(;
        N::Int = 10^2,
        N_prods::Int = 20, seed::Int = 0,
        prime :: Int128 = Int128(101))

    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end
    for _ in 1:N
        p_base = one(PolynomialModP, prime)
        for _ in 1:N_prods
            p = rand(PolynomialModP, prime)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_mod_p - PASSED")
end

"""
Test derivative of polynomials (as well as product).
"""
function prod_derivative_test_poly_mod_p(;
        N::Int = 10^2,
        seed::Int = 0,
        prime :: Int128 = Int128(101))
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_poly_mod_p - PASSED")
end


"""
Test division of polynomials modulo p.
"""
function division_test_poly_mod_p(;
        prime::Int128 = Int128(101),
        N::Int = 10^4,
        seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        p2 == 0 && continue
        p_prod = p1*p2
        q, r = zero(PolynomialModP, prime), zero(PolynomialModP, prime)
        try
            q, r = divide(p_prod, p2)
            (q, r) == (nothing,nothing) && continue
        catch e
            if typeof(e) == DivideError
                @assert p2 == 0
            else
                throw(e)
            end
        end
        @assert iszero(q*p2+r - p_prod)
    end
    println("division_test_poly_mod_p - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly_mod_p(;
        prime::Int128=Int128(3),
        N::Int = 10^3,
        seed::Int = 0)

    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        g, s, t = extended_euclid_alg(p1, p2)
        @assert s*p1 + t*p2 - g == 0
    end
    println("\next_euclid_test_poly_mod_p - PASSED")
end


"""
Super sparse polynomials

"""
function pow_mod_poly_mod_p(;
        prime::Int128 = Int128(19),
        N::Int = 30,
        seed::Int = 0)

    Random.seed!(seed)
    print("pow_mod_poly_mod_p\t")
    for k in 1:N
        f = rand(PolynomialModP, prime; prob_term = 0.01, mean_degree = Float64(120))
        @assert leading(f^k).coeff == (leading(f).coeff)^k
        print(".")
    end
    println("\tPASSED")
end


function chinese_remainder_multiplication_test(;N :: Int = 10^2, seed :: Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        f = rand(PolynomialSparse128)
        g = rand(PolynomialSparse128)
        @assert f * g == multiply(f, g)
    end
end
