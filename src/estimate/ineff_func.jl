function ineff(φn::Float64, φn_1::Float64, yt::Array{Float64, 1}, perror::Array{Float64}, H::Array{Float64},T::Int64)
    N = size(perror, 2)
    w=Array{Float64,1}
    # Inefficiency function during tempering steps
    if T==1
        for i=1:N
            w[i] = (φn/φn_1 )^(length(yt)/2)*exp(-1/2*perror(:,i)'*(φn-φn_1 )*inv(H)*perror(:,i))
        end
    #Inefficiency function during initialization (note φn_1 is never called)
    else
        for i=1:N
            w[i] = (φn/(2*pi))^(length(yt)/2)*det(H)^(-1/2)*exp((-1/2)*perror(:,i)'*φn*inv(H)*perror(:,i))
        end
    end
    W = w/mean(w)
    return sum(W.^2)/N
end