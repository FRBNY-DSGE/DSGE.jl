function mutation_RWMH1(m::AbstractModel, yt::Array{Float64,1},s_init, eps_init)

    #yt is a 1x3 vector of all three y at given time t
    #n_part = get_setting(m,:n_particles)
    N_MH=1
    #create covariance matrix of measurement error
    #H= diagm([m[:e_y].value,m[:e_π].value,m[:e_R].value])
    φ,R,C=solve(m)
    measure = measurement(m,φ,R,C)
    H=measure[:EE]
    varStateEq=measure[:MM]*measure[:QQ]*measure[:MM]'
    #data are wide so rows are vars an cols are time periods
    c=get_setting(m,:c)
    U,E,V=svd(varStateEq)
    cov_mat=U*diagm(sqrt(E))
    acpt=0
    for i=1:N_MH
        ##figure out what c should be
        eps_new=eps_init + c*cov_mat*randn(3)
        s_new_fore = φ*s_init+R*varStateEq*eps_new
        s_init_fore = φ*s_init+R*varStateEq*eps_init
        perror_new = yt-measure[:ZZ]*s_new_fore-measure[:DD]
        perror_init = yt-measure[:ZZ]*s_init_fore - measure[:DD]
        
        post_new = log(pdf(MvNormal(zeros(length(yt)),H),perror_new)*pdf(MvNormal(zeros(length(eps_new)),eye(length(eps_new),length(eps_new))),eps_new))
        post_init = log(pdf(MvNormal(zeros(length(yt)),H),perror_init)*pdf(MvNormal(zeros(length(eps_init)),eye(length(eps_init),length(eps_init))),eps_init))
       
        α = exp(post_new - post_init)
        
        if rand()<α 
            #accept
            s_init = s_new_fore
            eps_init = eps_new
            acpt = acpt+1
        #else reject and keep the old particle unchanged
        end
    end
    return s_init, eps_init, acpt 
end

