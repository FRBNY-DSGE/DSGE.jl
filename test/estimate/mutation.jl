using DSGE,DataFrames,HDF5

# Set model
m=AnSchorfheide(testing=true)
m<=Setting(:date_forecast_start,quartertodate("2015-Q4"))
m<=Setting(:tpf_N_MH,3)
# Set seeding 
srand(1234)
# Set path
path=dirname(@__FILE__)
# Compute system and store parameters
sys=compute_system(m)
Φ = sys.transition.TTT
A = sys.measurement.DD
B = sys.measurement.ZZ
R = sys.transition.RRR
H = sys.measurement.EE
varStateEq = sys.measurement.MM*sys.measurement.QQ*sys.measurement.MM'
U,E,V=svd(varStateEq)
cov_mat = U*diagm(sqrt(E))
sqrtS2=R*cov_mat'
N_MH = get_setting(m,:tpf_N_MH)

h5open("$path/../reference/matricesForMutation.h5","w") do file
    write(file, "phi", sys.transition.TTT)
    write(file, "D", sys.measurement.DD)
    write(file, "R", sys.transition.RRR)
    write(file, "H", sys.measurement.EE)
    write(file, "C", sys.transition.CCC)
    write(file, "M", sys.measurement.MM)
    write(file, "Q", sys.measurement.QQ)
    write(file, "Z", sys.measurement.ZZ)
    write(file, "N_MH", N_MH)
end
    

#Test that it compiles
s_part1, eps_part1, acpt = mutation(m,[50.2,8.3,7.6],[.8,.9,.6,.9,.11,5,7,10],[.2,.5,.7],A,B,R,Φ,H,sqrtS2,cov_mat,N_MH)

h5open("$path/../reference/mutation_RWMH1.h5","w") do file
    write(file, "s_part1", s_part1)
    write(file, "eps_part1",eps_part1)
end

data=h5open("$path/../reference/mutation_RWMH.h5","r") do file
    read(file,"data")
end

# Test function with one column of data.
s, eps, acpt = mutation(m, data[:,1], ones(8), zeros(3),A,B,R,Φ,H,sqrtS2,cov_mat,N_MH)  
c = h5open("$path/../reference/mutation_RWMH1.h5","r") do file
    read(file,"s_part1")
    read(file,"eps_part1")
end

# Test Kalman-certified reasonable case (should accept)
ind_s, ind_eps, ind_acpt = mutation(m, data[:,1], [-0.6119,-1.5262,-2.6471,0.9952,-0.8350,-0.4785,-0.3428,-1.6893], zeros(3),A,B,R,Φ,H,sqrtS2,cov_mat,N_MH)
@test ind_acpt == N_MH #very likely so should be equal to number of MH steps

#Test unreasonable case (should reject)
ind_s, ind_eps, ind_acpt = mutation(m, data[:,1], 3*[-0.6119,-1.5262,-2.6471,0.9952,-0.8350,-0.4785,-0.3428,-1.6893], zeros(3),A,B,R,Φ,H,sqrtS2,cov_mat,N_MH)
@test ind_acpt == 0

nothing 