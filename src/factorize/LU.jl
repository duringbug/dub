"""
struct LU{T<:Number}

这个方法对A进行LU分解,并且对A进行行变换,返回矩阵L与U

Arguments:
- `A`: 待分解的矩阵(n•n).
 
Returns:
- `A`: 分解后矩阵(n•n).经过了行变换
- `P`: 行变换矩阵
- `L`: L矩阵
- `U`: U矩阵

Examples:
```
julia> using dub

julia> lu = dub.LU([
           2//1 4//1
           4//1 5//1
       ])
dub.LU{Rational{Int64}}(Rational{Int64}[2//1 4//1; 4//1 5//1], 2, Rational{Int64}[1//1 0//1; 2//1 1//1], Rational{Int64}[2//1 4//1; 0//1 -3//1])

julia> println(lu.L * lu.U == lu.A)
true

julia> println(lu.L)
Rational{Int64}[1//1 0//1; 2//1 1//1]

julia> println(lu.U)
Rational{Int64}[2//1 4//1; 0//1 -3//1]
```


"""
struct LU{T<:Number}
    A::Matrix{T}
    n::Int
    L::Matrix{T}
    U::Matrix{T}
    P::Matrix{T}
    function LU{T}(A) where {T<:Number}
        n, A = check_square_matrix(A)
        P, A, L, U = lu_factorize(A, n, T)
        new{T}(A, n, L, U, P)
    end
    function LU(A::Matrix{T}) where {T<:Number}
        n, A = check_square_matrix(A)
        P, A, L, U = lu_factorize(A, n, T)
        new{T}(A, n, L, U, P)
    end

    function check_square_matrix(A::Matrix{T}) where {T<:Number}
        if size(A, 1) != size(A, 2)
            throw(ArgumentError("Matrix A must be square: $(size(A, 1)) != $(size(A, 2))"))
        end
        n = size(A, 1)
        return n, A
    end
    function lu_factorize(A::Matrix{T}, n, S) where {T<:Number}
        L = zeros(S, n, n)
        U = zeros(S, n, n)
        P = zeros(S, n, n)
        for i in 1:n
            P[i, i] = one(S)
        end
        A_i_tilde = A
        for i in 1:n
            A_i_tilde = P * (A - L * U)
            if A_i_tilde[i, i] == 0 && i != n
                for k in i+1:n
                    if A_i_tilde[k, i] != 0
                        A_i_tilde[i, :], A_i_tilde[k, :] = A_i_tilde[k, :], A_i_tilde[i, :]
                        P[i, :], P[k, :] = P[k, :], P[i, :]
                        break
                    end
                end
            end
            u_i = A_i_tilde[i, :]
            if A_i_tilde[i, i] == 0
                error("Matrix A_i_tilde is not full rank (singular).")
            end
            l_i = A_i_tilde[:, i] / A_i_tilde[i, i]
            U[i, :] = u_i
            L[:, i] = P * l_i
        end

        return P, P * A, P * L, U
    end
end