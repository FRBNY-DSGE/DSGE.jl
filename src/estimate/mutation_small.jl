
#function mutation(m::AbstractModel, system::System{Float64}, yt::Array{Float64,1}, s_init::Array{Float64,1}, ε_init::Array{Float64,1}, cov_s::Array{Float64,2}, nonmissing::Array{Bool,1})
#function mutation_problem(c::Float64, N_MH::Int64, deterministic::Bool, system::System{Float64}, yt::Array{Float64,1}, s_init::Array{Float64,1}, ε_init::Array{Float64,1}, cov_s::Array{Float64,2}, nonmissing::Array{Bool,1}, distCall::Bool)
#function mutation_problem2(c::Float64, N_MH::Int64, deterministic::Bool, system::System{Float64}, yt::Array{Float64,1}, s_init::AbstractArray{Float64,1}, ε_init::AbstractArray{Float64,1}, cov_s::Array{Float64,2}, nonmissing::Array{Bool,1}, distCall::Bool) 
function mutation_small(c::Float64, N_MH::Int64, deterministic::Bool, system::System{Float64}, yt::Array{Float64,1}, s::AbstractArray{Float64,2}, ε::AbstractArray{Float64,2},index::Int64, cov_s::Array{Float64,2}, nonmissing::Array{Bool,1}, distCall::Bool) 
   #------------------------------------------------------------------------
    # Setup
    #------------------------------------------------------------------------
    # Set path--used for testing
    path = dirname(@__FILE__)
    #s_init = s_i
    #ε_init = ε_i
    
    ##s_init = s[:,index]
    #@show s_init 
    ε_init = ε[:,index]
    # size(s)
    # size(ε)
    ##s_init = randn(54)
    #ε_init = randn(7)
    
    s_init = ones(size(s,1))
    for k=1:size(s,1)
        s_init[k]=s[k,index]
    end

    DD = system.measurement.DD[nonmissing]
    ZZ = system.measurement.ZZ[nonmissing,:]
    EE = system.measurement.EE[nonmissing,nonmissing]
    RRR = system.transition.RRR[:,nonmissing]
    TTT = system.transition.TTT
    QQ = system.measurement.QQ[nonmissing,nonmissing]
    sqrtS2 = RRR*Matrix(chol(nearestSPD(QQ)))'

    # Initialize ind_s and ind_ε
    ind_s = s_init
    ind_ε = ε_init
    
 
    #------------------------------------------------------------------------
    # Metropolis-Hastings Steps
    #------------------------------------------------------------------------
    for i=1:N_MH
        rand_mat = randn(size(QQ,1),1) 
        ε_new=ε_init + c*Matrix(chol(nearestSPD(cov_s)))'*rand_mat

        s_new_fore = TTT*s_init+sqrtS2*ε_new

        s_init_fore = TTT*s_init+sqrtS2*ε_init

        error_new  = yt - ZZ*s_new_fore - DD
        error_init = yt - ZZ*s_init_fore - DD

        # Calculate the top and bottom probabilities for the α ratio.
       if distCall
           post_new = log(pdf(MvNormal(zeros(length(yt)),EE),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
           post_init = log(pdf(MvNormal(zeros(length(yt)),EE),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])

        # α is the probability of accepting the new particle 
        # (post_new and post_init are in logs so subtract and take the exponential)
        α = exp(post_new - post_init)
       else
           α = .5
       end
        # Accept the particle with probability α
        if (deterministic) num = 0.5
        else num = rand() end
        
        if num<α 
            # Accept
            ind_s = s_new_fore
            ind_ε = ε_new
        else 
            # Reject and keep the old particle unchanged
            ind_s = s_init_fore
            ind_ε = ε_init
        end
        ε_init = ind_ε
    end   
    return ind_s, ind_ε
    gc()
end
