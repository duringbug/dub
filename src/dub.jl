module dub

function test()
    using .factorize.cholesky
    cholesky_result = cholesky([4 -1 -1; -1 4.25 4.25; 1 2.75 3.5])
    println(cholesky_result)
end

end # module dub