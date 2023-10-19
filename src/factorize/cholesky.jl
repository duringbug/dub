module function cholesky(A)
    n = size(A, 1)
    L = zeros(n, n)

    for j in 1:n
        L[j, j] = sqrt(A[j, j] - sum(L[j, 1:j-1]^2))
        for i in j+1:n
            L[i, j] = (A[i, j] - sum(L[i, 1:j-1] * L[j, 1:j-1])) / L[j, j]
        end
    end
    return L
end
end