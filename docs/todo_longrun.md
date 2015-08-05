# Long-Run Things to Do

## Efficiency
- Reimplement eqcond, transition, and measurement matrices using sparse matrices
- Add option for `eqcond` to take in Γ0, Γ1, etc. matrices as arguments so we can reassign to those matrices instead of reallocating `zeros(82, 82)` each time
- Same with other functions that also do this
- Replace csminwel, Hessian computation, etc. with external packages

## Organization
