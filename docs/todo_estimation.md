# Step 1

## Pearl
- confirm prior density calculated in Julia matches prior density calculated in Matlab
  - works now
  - note that we had X~IG(\alpha,\beta) but we actually want X^2~IG(\alpha,\beta)
- finish calculating likelihood/check objective function [hardest]

## Erica
- csminwel (should be completely self-contained)

# Step 2: Compute the inverse hessian

## Pearl
- if we have a covariance matrix saved, calculated SVD of that
- else, use one of two things:
  - Julia port of hessizero (preferred) [may be difficult, but hopefully not]
  - existing function, such as from ForwardDiff.jl package (as time allows,
    i.e. once estimation is complete)

## Erica
- compute inverse hessian via SVD [hopefully straightforward. make sure
  Matlab/Julia SVD routines are analogous]

# Step 3: Actual MH algorithm

- initialize by finding a valid parameter vector
- draw from the proposal distribution [straightforward]
- solve the model [done]
- compute the likelihood again [done]
- choose to accept/reject [straightforward]
- save objects of interest [coding]
- calculate covariance matrix [straightforward]
