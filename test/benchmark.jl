using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")


# Benchmark CRT multiplication
function chinese_remainder_multiplication(;N :: Int = 10^4, seed :: Int = 0)
    println("Old Multiplication")
    @time benchmark(*, N, seed)
    println("CRT Multiplication")
    @time benchmark(multiply, N, seed)
end

function benchmark(fn, N :: Int, seed :: Int)
    Random.seed!(seed)
    for _ in 1:N
        f = rand(PolynomialSparse128)
        g = rand(PolynomialSparse128)
        @assert f * g == fn(f, g)
    end
end


chinese_remainder_multiplication()
