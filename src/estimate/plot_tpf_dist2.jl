using DSGE, HDF5, DataFrames
using QuantEcon: solve_discrete_lyapunov
using Plots
include("tpf_error.jl")

path = dirname(@__FILE__)

# error = zeros(10)

# m = SmetsWouters("ss1", testing=true)
# #m = AnSchorfheide(testing=true)

# for i=1:10
#     error[i] = tpf_error(m)
#     if typeof(m)==SmetsWouters{Float64}
#         h5open("$path/../../test/reference/error_block11_sw.h5","w") do file
#             write(file, "error", error)
#         end
#     else
#         h5open("$path/../../test/reference/error_block10_as.h5","w") do file
#             write(file, "error", error)
#         end
#     end
# end


errors = zeros(0)
for i=1:10
    error = h5read("$path/../../test/reference/error_block$(i)_as.h5","error")
    errors = vcat(errors,error)
end
# for i=[12,13,14,15,17,18,19]
#     error = h5read("$path/../../test/reference/error_block$(i)_sw.h5","error")
#     errors = vcat(errors,error)
# end
errors = errors[errors.!=0]
#this line fixes to match paper because i subtracted in the wrong order
errors = (-1).*errors
@show length(errors)
plotly()
histogram(errors)
#plot!(title = "a")
gui()
