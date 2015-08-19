import DSGE
# import Hobo

# Example usage:
function rosenbrock(x::Vector)
    (1-x[1])^2.0 + 105*(x[2]-x[1]^2.0)^4.0
end

# this is the actual gradient of the Rosenbrock function
function rosenbrock_grad(x::Array)
    dr = similar(x)
    dr[1] = -2*(1 - x[1]) - 8*105*x[1]*(x[2]-x[1]^2.0)^3.0
    dr[2] = 4*105*(x[2]-x[1]^2.0)^3.0
    badg = false
    return dr, badg
end

# A really bad guess
x_init = [10.0, -9.0]

# println("Testing using Hobo.csminwel...")
# res_real_grad  = Hobo.csminwel(rosenbrock, rosenbrock_grad, [10.0, -9.0])
# res_numeric_grad  = Hobo.csminwel(rosenbrock, [10.0, -9.0])

res_real_grad  = DSGE.csminwel(rosenbrock, rosenbrock_grad, [10.0, -9.0])
res_numeric_grad  = DSGE.csminwel(rosenbrock, [10.0, -9.0])

None
