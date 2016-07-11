function mutation_RWMH(p0, l0, post0, tune, i, f)
# RWMH for mutation step
# INPUT
# tune
# p0 = para(j,:)'
# l0 = loglh(j)
# post0 = logpost(j)

# OUTPUT
# ind_para
# ind_loglh
# ind_post
# ind_acpt


# RW proposal
px = p0 + tune.c*chol(tune.R)'*randn(tune.npara,1)
[postx, lx] = f(px, tune.phi(i))

# Previous posterior needs to be updated (due to tempering)
post0 = post0+(tune.phi(i)-tune.phi(i-1))*l0

# Accept/Reject
alp = exp(postx - post0) # this is RW, so q is canceled out
if rand < alp # accept
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
