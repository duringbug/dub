include("factorize/LU.jl")
lu = LU([
    2//1 4//1
    4//1 5//1
])
println(lu.L)
println(lu.U)
println(lu.L * lu.U == lu.A)