using DSGE,DataFrames,HDF5

m=AnSchorfheide(testing=true)
m<=Setting(:date_forecast_start,quartertodate("2015-Q4"))

s_part1, eps_part1, acpt = mutation_RWMH1(m,[50.2,8.3,7.6],[.8,.9,.6,.9,.11,5,7,10],[.2,.5,.7])

path=dirname(@__FILE__)
h5open("$path/../reference/mutation_RWMH1.h5","w") do file
    write(file, "s_part1", s_part1)
    write(file, "eps_part1",eps_part1)
end

c = h5open("$path/../reference/mutation_RWMH1.h5","r") do file
    read(file,"s_part1")
    read(file,"eps_part1")
end
println(c)

#path=dirname(@__FILE__)
#acc_log1h= h5read("$path/../reference/mutation_RWMH.h5","acc_log1h")

#=
df,system,s0,eps0 = jldopen($"path/../reference/forecast_args.jld","r") do file
    read(file,"df"),read(file,"system"),read(file,"s0",read(file,"eps0"
end

#read expected output
exp 
=#