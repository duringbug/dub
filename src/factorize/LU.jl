```
对矩阵A进行LU分解,采用的方法为秩一分解法
```
struct LU{T<:Number}
    A::Matrix{T}
    n::Int
    L::Matrix{T}
    U::Matrix{T}
    function LU{T}(A::Matrix{T}) where {T<:Number}
        n, A = check_square_matrix(A)
        L, U = lu_factorize(A, n)
        new{T}(A, n, L, U)
    end
    function LU(A::Matrix{T}) where {T<:Number}
        n, A = check_square_matrix(A)
        L, U = lu_factorize(A, n)
        new{T}(A, n, L, U)
    end

    function check_square_matrix(A::Matrix{T}) where {T<:Number}
        if size(A, 1) != size(A, 2)
            throw(ArgumentError("Matrix A must be square: $(size(A, 1)) != $(size(A, 2))"))
        end
        n = size(A, 1)
        return n, A
    end
    function lu_factorize(A::Matrix{T}, n) where {T<:Number}
        L = zeros(T, n, n)
        U = zeros(T, n, n)
        A_i_tilde = A
        for i in 1:n
            A_i_tilde = A - L * U
            u_i = A_i_tilde[i, :]
            l_i = A_i_tilde[:, i] / A_i_tilde[i, i]
            U[i, :] = u_i
            L[:, i] = l_i
        end
        return L, U
    end
end