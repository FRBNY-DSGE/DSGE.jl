using DSGE,DataFrames,HDF5

m=AnSchorfheide(testing=true)
m<=Setting(:date_forecast_start,quartertodate("2015-Q4"))
srand(1234)
path=dirname(@__FILE__)
sys=compute_system(m)
h5open("$path/../reference/matricesForMutation.h5","w") do file
    write(file, "phi", sys.transition.TTT)
    write(file, "D", sys.measurement.DD)
    write(file, "R", sys.transition.RRR)
    write(file, "H", sys.measurement.EE)
    write(file, "C", sys.transition.CCC)
    write(file, "M", sys.measurement.MM)
    write(file, "Q", sys.measurement.QQ)
    write(file, "Z", sys.measurement.ZZ)
end
    

#test it. it works!
s_part1, eps_part1, acpt = mutation(m,[50.2,8.3,7.6],[.8,.9,.6,.9,.11,5,7,10],[.2,.5,.7])


h5open("$path/../reference/mutation_RWMH1.h5","w") do file
    write(file, "s_part1", s_part1)
    write(file, "eps_part1",eps_part1)
end

data=h5open("$path/../reference/mutation_RWMH.h5","r") do file
    read(file,"data")
end



#test it with one column of data. it works!
#s,eps,acpt=mutation(m,data[:,1],[.8,.9,.6,.9,.11,5,7,10],[.2,.5,.7])
 s, eps, acpt = mutation(m, data[:,1],ones(8), zeros(3))  
print(s)
c = h5open("$path/../reference/mutation_RWMH1.h5","r") do file
    read(file,"s_part1")
    read(file,"eps_part1")
end

#path=dirname(@__FILE__)
#acc_log1h= h5read("$path/../reference/mutation_RWMH.h5","acc_log1h")

#=
df,system,s0,eps0 = jldopen($"path/../reference/forecast_args.jld","r") do file
    read(file,"df"),read(file,"system"),read(file,"s0",read(file,"eps0"
end

#read expected output
exp 
=#