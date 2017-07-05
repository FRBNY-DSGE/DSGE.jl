"""
### Arguments:

- 'm': model object
- 'data': a matrix or data frame containing time series of the observation

### Output:
- gotta figure that out.
"""
using Distributions

function mutate(m::AbstractModel, data::Matrix, ϵ_0::Vector{AbstractFloat}, s0::Vector{AbstractFloat}, t::Int64, verbose::Symbol=:low)
    
    ind_acpt{Int64} = 0
    N_mh = get_setting(m, :n_mh_pf)
    c = get_setting(m, :c)

    Φ, RRR, CCC = solve(m)
    ### m_out = measurement(m, Φ, RRR, CCC)
    m_sys = compute_system(m)
    HH = m_sys.measurement.EE + m_sys.measurement.MM*m_sys.measurement.QQ*m_sys.measurement.MM'
    ### HH = m_out[:EE]+m_out[:MM]*m_out[:QQ]*m_out[:MM]'
    
    sqrtS2 = RRR*svd(m_out[:QQ])[1]'

    for nn = 1:N_mh

        ϵ_x = ϵ_0 + c*svd(RRR)[1]'*rand(size(RRR,1),1)

        sx_fore = Φ*s0 + sqrtS2*ϵ_x
        s0_fore = Φ*s0 + sqrtS2*ϵ_0

        perrorx = yy'-A-B*sx_fore
        perror0 = yy'-A-B*s0_fore

        postx = log(pdf(MvNormal(zeros(1,size(yy,2)),HH),perrorx')*pdf(MvNormal(ϵ_x')))
        post0 = log(pdf(MvNormal(zeros(1,size(yy,2)),HH),perror0')*pdf(MvNormal(ϵ_0'))))

        alp = exp(postx-post0)
        if rand < alp:
            ind_s = sx_fore
            ind_ϵ = ϵ_x
            ind_acpt = ind_acpt+1
        else
            ind_s = s0_fore
            ind_ϵ = ϵ_x
            #ind_acpt = 0
        end

        ϵ_0 = ind_ϵ
    end
    ind_acpt = ind_acpt/N_mh
    return [ind_s, ind_ϵ, ind_acpt]
end