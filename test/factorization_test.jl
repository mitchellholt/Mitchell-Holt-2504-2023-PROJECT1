#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################


"""
Test factorization of polynomials.
"""
function factor_test_poly(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialDense)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly - PASSED")
end

function factor_test_poly_sparse(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_sparse - PASSED")
end

function factor_test_poly_sparse_128(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse128)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_sparse_128 - PASSED")
end

function factor_test_poly_mod_p(;
        N::Int = 10,
        seed::Int = 30,
        primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialModP, Int128(prime))
            factorization = factor(p)
            pr = expand_factorization(factorization)
            @assert p - pr == 0
        end
    end

    println("\nfactor_test_poly_mod_p - PASSED")
end


function factor_test(;
        N :: Int = 10,
        prime :: Int128 = Int128(19),
        deg :: Int = 3,
        seed :: Int = 0)

    Random.seed!(seed)
    for _ in 1:10
        f = rand(PolynomialModP, prime, degree = deg)
        g = rand(PolynomialModP, prime, degree = deg)
        factorization = factor(f * g)
        pr = expand_factorization(factorization)
        @assert f * g - pr == 0
        print(".")
    end
    println()
end
