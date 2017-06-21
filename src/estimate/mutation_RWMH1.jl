functiomn mutation_RWMH1(m::AbstractModel, yt::Matrix{Float64}, R, c, S0)

    #yt is a 1xn vector of all y for a given time
    n_part = get_setting(m,:n_particles)
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

    for i=1:N_MH
        ##figure out what c should be
        eps_new=eps_init + c*cov_mat*randn(3,1)
        s_new_fore = φ*s_init+R*varStateEq*eps_new
        s_init_fore = φ*s_init+R*varStateEq*eps_init
        perror_new = data-measure[:ZZ]*s_new-measure[:DD]
        perror_init = data-measure[:ZZ]*s_init - measure[:DD]
        
        post_new = log(pdf(MvNormal(zeros(length(yt)),H), perror_new)*pdf(MvNormal(eps_new))
        post_init = log(pdf(MvNormal(zeros(length(yt)),H), perror_init)*pdf(MvNormal(eps_init))
        α = exp(pos_new - post_init)
        if rand()<α 
            #accept
            s_init = s_new
            eps_init = eps_new
            acpt = acpt+1
        #else reject and keep the old particle unchanged
        end
    end
    return s_init, eps_init, acpt 
end