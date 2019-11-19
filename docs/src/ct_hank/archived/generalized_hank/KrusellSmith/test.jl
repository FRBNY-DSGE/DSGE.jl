using ForwardDiff
ForwardDiff.jacobian(x -> ForwardDiff.jacobian(cumprod, x), [1,2,3])


