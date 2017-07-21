
function pmap_test1(x::Int64)
    ans = factorial(round(x/10))
    sleep(.5)
    return ans,2,3
end

#Pmap faster than parallel than no parallel
function pmap_test2(A,B,C,D)
    ans1=0
    ans2=0
    for i = 1:10
        ans1 = A*B
        ans1 = ans1'
        ans1 = (ans1*randn(500,1))'*randn(500,1)
        ans2 = C*D
        ans2 = ans2'
        ans2 = (ans2*randn(1000,1))'*randn(700,1)
     end
     return ans1[1],ans2[1], 3
end

function pmap_test3(A,B,C,D)
    ans1 = A*B
    ans1 = ans1'
    ans1 = (ans1*randn(500,1))'*randn(500,1)
    
    ans2 = C*D
    ans2 = ans2'
    ans2 = (ans2*randn(1000,1))'*randn(700,1)
    return ans1[1],ans2[1], 3
end

function pmap_test4(yt::Array{Float64,1}, system, s_init::Array{Float64,1}, ε_init::Array{Float64,1}, cov_s::Array{Float64,2}, nonmissing::Array{Bool,1})
    #path = dirname(@__FILE__)
    D = system.measurement.DD[nonmissing]
    Z = system.measurement.ZZ[nonmissing,:]
    E = system.measurement.EE[nonmissing,nonmissing]
    R = system.transition.RRR[:,nonmissing]
    T = system.transition.TTT
    Q = system.measurement.QQ[nonmissing,nonmissing]
    sqrtS2 = R*Matrix(chol(nearestSPD(Q)))'
    acpt=0
    ind_s = s_init
    ind_ε=ε_init
    for i=1:2
        rand_mat = randn(size(Q,1),1)
        
        ε_new = ε_init + .2*Matrix(chol(nearestSPD(cov_s)))'*rand_mat
        
        s_new = T*s_init+sqrtS2*ε_new
        s_init = T*s_init+sqrtS2*ε_init
        
        error_new = yt - Z*s_new - D
        error_init = yt - Z*s_init - D
        
        post_new = log(pdf(MvNormal(zeros(length(yt)),E),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(length(yt)),E),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])
        α = exp(post_new - post_init)

        if rand()<α
            ind_s = s_new
            ind_ε=ε_new
            acpt += 1
        else
            ind_s = s_init
            ind_ε=ε_init
        end
        ε_init = ind_ε
    end
    acpt/=2
    return ind_s, ind_ε, acpt
end

using DSGE, Optim, Distributions