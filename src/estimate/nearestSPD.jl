"""
```
nearestSPD(A::Matrix)
```
Calculates the nearest (in Frobenius norm) symmetric positive semidefinite matrix to A. 

From Higham: The nearest symmetric positive semidefinite matrix in the
Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
 where H is the symmetric polar factor of B=(A + A')/2.

http://www.sciencedirect.com/science/article/pii/0024379588902236

### Inputs

- `A::Matrix`: A square matrix

### Outputs
 
- `Ahat::Matrix`: The symmetric, square positive semidefinite matrix closest to A
"""
function nearestSPD(A::Matrix)
     
    if size(A,1) != size(A,2) 
        error("A must be a square matrix")
    end
    # symmetrize A into B
    B = (A + A')/2;
    
    # Compute the symmetric polar factor of B. Call it H.
    # Clearly H is itself SPD.
    U, Σ, V = svd(B);
    H = V*diagm(Σ)*V';
    
    Ahat = (B+H)/2;
    # ensure symmetry
    Ahat = (Ahat + Ahat')/2;
    
    # test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    k = 1 
    while minimum(eigvals(Ahat)) < 0 
        mineig = minimum(eigvals(Ahat))
        Ahat = Ahat + (-mineig*k^2 + eps(mineig))*eye(size(A,1),size(A,2))
        k = k + 1
    end
    
    return Ahat
end