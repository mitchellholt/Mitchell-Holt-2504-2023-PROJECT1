using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


function show(io :: IO, factors :: Vector{Tuple{PolynomialSparse, Int}})
    for (factor, m) in factors
        print(io, "($factor)", m == 1 ? "" : "^$(int_to_superscript(m))")
    end
end

x = x_poly(PolynomialSparse)
Random.seed!(0)

# Addition of polynomials
println("Basic arithemtic operations on random polynomials")
p1, p2 = rand(PolynomialSparse), rand(PolynomialSparse)
@show p1
@show p2
@show p1 + p2

# Product of polynomials
@show p1 * p2
println("\nProduct of polynomials")
@show 25 * (x^2 + 1) * (x^2 + (-1))
@show (x + 1)^9

# extended euclidean algorithm
println("\nEuclidean division of two polynomials in Z/19[x]")
prime = 19
@show prime
q1, q2 = (x^2 + 1)*(x^2 + -1), (x + 1)*(x + -1)
@show q1
@show q2
(gcdiv, s, t) = extended_euclid_alg(q1, q2, prime)
println("gcd(q1, q2, prime) = $gcdiv")
println("Construction of the gcd using Bezierre Coefficients s*q1 + t*q2")
@show s
@show t
@show mod(q1*s + q2*t, prime) == mod(gcdiv, prime)

# division
println("\nDivision")
f1 = rand(PolynomialSparse; degree = 4, mean_degree = 2.0, max_coeff = 5)
f2 = rand(PolynomialSparse; degree = 3, mean_degree = 1.5, max_coeff = 5)
f3 = rand(PolynomialSparse; degree = 2, mean_degree = 1.0, max_coeff = 5)
g1, g2 = f1 * f2 + f3, f2
@show g1
@show g2
println("Euclidean division for g1 / g2:")
primes = [7, 11, 13, 17, 19, 23, 29]
for p in primes
    dividend, remainder = (divide(g1, g2))(p)
    println("\tin Z/$p[x]:\tg1 = ($dividend)($g2) + $remainder  (mod $p)")
end

# factorization
println("\nFactorization in Z/p[x]")
println("Let f = $f1 * $f2")
f = f1 * f2
@show f
for p in primes[3:end]
    factorization = factor(f, p)
    println("\tin Z/$p[x]:\tf = $factorization")
    println("\tReconstruction:\t$(mod(expand_factorization(factorization), p))")
end
