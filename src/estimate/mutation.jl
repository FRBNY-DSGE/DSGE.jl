"""
```
function mutation(m::AbstractModel, yt::Array{Float64,1},s_init::Array{Float64,1}, ε_init::Array{Float64,1}, A::Array, B::Array, R::Array, Φ::Array, H::Array, sqrtS2::Array, cov_mat::Array,N_MH::Int64)
```
Runs random-walk Metropolis Hastings for ONE particle. The caller should loop through all particles, calling this method on each. 
m: The model. Used to get the setting tpf_c (which is VERY important it turns out...) 
yt: A 1 x num_measurement vector at time t for all observed y (GDP, PCE, FFR, etc.)
s_init: The starting state before mutation.
ε_init: The starting epsilon (state error) before mutation.
"""

function mutation(m::AbstractModel, yt::Array{Float64,1},s_init::Array{Float64,1},ε_init::Array{Float64,1},cov_s::Array,rand_mat::Array, testing::Int64)
    #------------------------------------------------------------------------
    # Setup
    #------------------------------------------------------------------------
    # Set path--used only for testing
    path = dirname(@__FILE__)

    A = get_setting(m, :DD)
    B = get_setting(m, :ZZ)
    R = get_setting(m, :RRR)
    Φ = get_setting(m, :TTT)
    H = get_setting(m, :EE)
    sqrtS2 = get_setting(m, :tpf_sqrtS2)
    N_MH = get_setting(m, :tpf_N_MH)

    #Initialize ind_s and ind_ε
    ind_s=s_init
    ind_ε=ε_init
    #Get the meta-paramater tpf_c from model setting. This updates dynamically to tune TPF over time (which is very important as we all know).
    c=get_setting(m,:tpf_c)
    # Initialize acceptance counter to zero
    acpt=0
    
    #------------------------------------------------------------------------
    # Metropolis-Hastings Steps
    #------------------------------------------------------------------------
    for i=1:N_MH
        # Generate new draw of ε from a N(ε_init, c²R) distribution where R is cov_s (defined in tpf.jl) and c is the tuning parameter
        if testing==0 rand_mat=randn(3,1) end
        ε_new=ε_init + c*Matrix(chol(nearestSPD(cov_s)))'*rand_mat

        # Use the state equation to calculate the corresponding state from that ε 
        s_new_fore = Φ*s_init+sqrtS2*ε_new

        # Use the state equation to calculate the state corresponding to ε_init
        s_init_fore = Φ*s_init+sqrtS2*ε_init

       # Calculate difference between data and expected y from measurement equation and calculated states from above for both the new draw of ε_new and the old ε_init (we do this to calculate probabilities below. Since the error is still drawn from a Normal and everything is still linear, we know that y will also be normal. See equation of multivariate normal for how error_new and error_init enter into the pdf).
        error_new = yt-B*s_new_fore-A #wrong
        #s_new_fore is wrong i think because eps_new is wrong
        
        error_init = yt-B*s_init_fore-A
        # Calculate the top and bottom probabilities for the α ratio.
        post_new = log(pdf(MvNormal(zeros(length(yt)),H),error_new)[1]*pdf(MvNormal(zeros(length(ε_new)),eye(length(ε_new),length(ε_new))),ε_new)[1])
        post_init = log(pdf(MvNormal(zeros(length(yt)),H),error_init)[1]*pdf(MvNormal(zeros(length(ε_init)),eye(length(ε_init),length(ε_init))),ε_init)[1])

        # α represents the probability of accepting the new particle (post_new and post_init are in logs so subtract and take the exponential)
        α = exp(post_new - post_init)

        # Accept the particle with probability α
        if (testing==0) num = rand()
        else num=0.5 end
        if num<α 
            # Accept
            ind_s = s_new_fore
            ind_ε = ε_new
            acpt = acpt+1
        else 
            # Reject and keep the old particle unchanged
            ind_s = s_init_fore
            ind_ε = ε_init
            #THIS IS DIFFERENT THAN ORIGINAL MATLAB BUT WE THINK MATLAB WRONG
            #I HAVE RE-UNCOMMENTED THIS SO THAT C DOESNT GET HUGE AND EVERYTHING DOESN'T BREAK WITH THE INDEX OUT OF BOUNDS ERROR
            if testing==0 acpt = 0 end
        end
        ε_init = ind_ε
    end
    acpt /= N_MH
    return ind_s, ind_ε, acpt 
end

