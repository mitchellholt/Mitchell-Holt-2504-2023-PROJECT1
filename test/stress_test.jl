using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")


function factor_test(;
        N :: Int = 10,
        prime :: Int128 = Int128(19),
        deg :: Int = 3,
        seed :: Int = 0)

    Random.seed!(seed)
    while true
        f = rand(PolynomialModP, prime, degree = deg)
        g = rand(PolynomialModP, prime, degree = deg)
        factorization = factor(f * g)
        pr = expand_factorization(factorization)
        @assert f * g - pr == 0
        print(".")
    end
    println()
end


factor_test()
