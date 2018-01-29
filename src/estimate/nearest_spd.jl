"""
```
nearest_spd(A::Matrix)
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
function nearest_spd(A::Matrix)

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
    # Ensure symmetry
    Ahat = (Ahat + Ahat')/2;

    # Test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    k = 0
    chol_success = false

    eig_matrix = eigvals(Ahat)
    while !chol_success
        try
            chol(Ahat)
            chol_success = true
        end
        k = k + 1
        if !chol_success
            # mineig_factor is used in the event that the matrix has very small, slightly
            # negative eigenvalues, then the identity scaled by the absolute value of the eigenvalue must
            # be added so as to make the matrix "more" PD and not less.
            # If the minimum eigenvalue is 0 then mineig_factor will permanently remain at 0
            # and the Ahat matrix will cease to change in the iteration, so the ternary
            # condition sets mineig_factor to be the machine epsilon for Float64 numbers
            # so that the algorithm will continue to perturb Ahat in the right direction and not
            # get stuck.

            eig_list = eigvals(Ahat)
            mineig = eig_list[1]
            mineig_factor = mineig != 0 ? abs(mineig) : eps(Float64)
            Ahat = Ahat + (mineig_factor*k^2 + eps(Float64))*eye(size(A)...)
            eig_matrix = hcat(eig_matrix, eigvals(Ahat))
        end
    end

    return Ahat
end
