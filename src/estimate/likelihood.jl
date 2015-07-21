# This is a dsge likelihood function that can handle 2-part estimation where
# there is a model switch.
# If there is no model switch, then we filter over the main sample all at once.
function likelihood{T<:FloatingPoint}(model::AbstractModel, YY::Array{T, 2})
    
    spec = model.spec_vars



    # TODO: Find a way to return these matrices from measurement equation or something so `solve` isn't called twice
    ## step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
    ## transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
    TTT, CCC, RRR = solve(model::AbstractModel)


    
    ## step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
    ## where u_t = eta_t+MM* eps_t with var(eta_t) = EE
    ## where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'

    # Get measurement equation matrices set up
    #try
        ZZ, DD, QQ, EE, MM = measurement(model)
    #catch 
        # Error thrown during gensys
    #    return -Inf
    #end

    HH = EE + MM*QQ*MM'
    VV = QQ*MM'
    VVall =  [[RRR*QQ*RRR' RRR*VV]; [VV'*RRR' HH]]

    # TODO: What's happening here? In the Matlab code, this was only done for period 2
    if any(CCC != 0)
        DD = DD + ZZ*((eye(spec["n_states_aug"]) - TTT)\CCC)
    end


    
    ## step 3: compute log-likelihood using Kalman filter - written by Iskander
    ##         note that Iskander's program assumes a transition equation written as:
    ##         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
    ##         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
    ##         and  VV2 = cov(eps2_t,u_u) = RRR*VV
    ##         define VVall as the joint variance of the two shocks VVall = var([eps2_tu_t])

    # TODO: Can we solve Lyapunov equation on matrices with anticipated shocks?
    # Solve lyapunov with normal period state matrices (i.e. period 2 matrices)
    A0 = zeros(size(TTT, 1), 1)
    P0 = dlyap(TTT, RRR*QQ*RRR')

    pyt, zend, Pend = kalcvf2NaN(YY', 1, zeros(spec["n_states_aug"], 1), TTT, DD, ZZ, VVall, A0, P0, 2)
    return pyt
end
