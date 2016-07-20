"""
'''
mutation_RWMH(p0, l0, post0, tune,i, data; rvec = [], rval = [], px = [], lx = [], postx = [])
'''

Execute one proposed move of the Metropolis-Hastings algorithm for a given parameter

### Arguments:
- 'p0': initial prior draw. para(j,:)
- 'l0': initial log-likelihood at the prior draw. loglh(j)
- 'post0': initial log-posterior at the prior draw. logpost(j)
- 'tune': the tune object, which is a struct/type object that contains information about tuning the SMC and the MH algorithms
- 'i': the index of the iteration of the SMC algorithm
- 'data': well-formed data as DataFrame

### Optional Arguments: 
- 'rvec': A matrix of horizontally concatenated random vectors ε as in Θ* = θ + ε for the proposed move in RWMH, for the purposes of testing that mutation_RWMH.
- 'rval': A vector of random values generated in MATLAB for the purposes of testing whether the proposed move exceeds the threshhold and is accepted (or rejected).
- 'px', 'lx', 'postx': Proposed new prior draw, log-likelihood, and log-posterior from executing the MATLAB code

### Output:

- 'ind_para': 
- 'ind_loglh': 
- 'ind_post': 
- 'ind_acpt':

"""
include("posterior.jl") #be sure to remove the call of f in line 20 since obj_func is 
#being replaced by posterior() in the new code base

function mutation_RWMH(p0, l0, post0, tune, i, data; rvec = [], rval = [], px = [], lx = [], postx = [])


# If testing, set the random seeds at fixed numbers
rvec = randn(tune.npara,1)
rval = rand()

if s.testing
    #srand(s.rng, 654)
    rvec = 
    rval =
end

#RW Proposal
px = p0 + tune.c*chol(tune.R)'*rvec
postx,lx = posterior(s,data)
    
# Previous posterior needs to be updated (due to tempering)
post0 = post0+(tune.phi(i)-tune.phi(i-1))*l0
# Accept/Reject
alp = exp(postx - post0) # this is RW, so q is canceled out
    
if rval .< alp # accept
    ind_para   = px
    ind_loglh  = lx
    ind_post   = postx
    ind_acpt   = 1
else
    ind_para   = p0
    ind_loglh  = l0
    ind_post   = post0
    ind_acpt   = 0
end

# # outside of function
# parasim(i,j,:) = ind_para
# loglh(j)       = ind_loglh
# temp_acpt(j,1) = ind_acpt
return ind_para, ind_loglh, ind_post, ind_acpt
end
