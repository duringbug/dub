function test()
    # lu = LU{Float32}([
    #     2//1 4//1
    #     4//1 5//1
    # ])
    lu = LU{Float64}([
        5//1 3//1 -1//1 3//1 5//1
        -5//1 -3//1 4//1 -4//1 2//1
        0//1 1//1 1//1 -2//1 1//1
        0//1 1//1 1//1 0//1 6//1
        0//1 1//1 1//1 3//1 3//1
    ])
    println(lu.L)
    println(lu.U)
    println(lu.L * lu.U)
    println(lu.L * lu.U == lu.A)
    qr = QR{Float64}([
        1.0 1.0 2.9
        2.0 0.0 3.4
        2.0 1.0 5.4
        3.5 1.2 3.1
    ])
    println(qr.Q)
    println(qr.R)
    println(qr.Q * qr.R)
    println(qr.Q * qr.R ≈ qr.A)
end