"""
Householder: Householder 镜像变换

GramSchmidt: Gram-Schmidt 正交化

Givens: Givens 旋转变换

"""
@enum QRType begin
    Householder = 1
    GramSchmidt = 2
    Givens = 3
end
"""
function QR(A, qr_type::QRType)

这个方法对A进行QR分解,返回矩阵Q与R

Arguments:
- `A`: 待分解的矩阵(M•N).
- `qr_type`: 分解的方法: @enum QRType Householder GramSchmidt Givens.
 
Returns:
QR分解的Q(M•M)与R(M•N)矩阵

Examples:
```
QR([4 -1 -1; -1 17//4 17//4; 1 11//4 7//2], Givens)
```
"""
function QR(A, qr_type::QRType)
    if Int(qr_type) == 3
        for j in 1:size(A, 2)
            x = view(A, :, j)
            println(x)
        end
    end
end