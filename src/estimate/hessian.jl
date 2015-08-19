# Compute Hessian of posterior function evaluated at x (vector)
# if noisy, display error messages, results, etc.
# 11/12/01 translated by Marco DelNegro in matlab from Frank Schorfheide's program in gauss
function hessizero!{T<:FloatingPoint}(x::Vector{T}, model::AbstractDSGEModel, YY::Matrix{T}; noisy::Bool = false)

    update!(model, x)

    ## index of free parameters
    para_free = [!α.fixed for α in model.Θ]
    fpara_free = find(para_free)
    nfree = length(fpara_free)

    npara = length(x)
    ndx = 6
    dx =  exp(-(6:2:(6+(ndx-1)*2))')
    hessian = zeros(npara, npara)
    #gradx = zeros(ndx, 1)
    #grady = zeros(ndx, 1)
    #gradxy = zeros(ndx, 1)
    hessdiag = zeros(ndx, 1)
    dxscale = ones(npara, 1)



    # Compute Diagonal elements first
    for seli = fpara_free'
        if noisy
            println("Hessian element: ($seli, $seli)")
        end

        for k = 1:ndx
            paradx = copy(x)
            parady = copy(x)
            paradx[seli] = paradx[seli] + dx[k]*dxscale[seli]
            parady[seli] = parady[seli] - dx[k]*dxscale[seli]
            #paradxdy = copy(paradx)
            #paradxdy[seli] = paradxdy[seli] - dx[k]*dxscale[seli]

            fx = posterior!(x, model, YY)
            fdx = posterior!(paradx, model, YY)
            fdy = posterior!(parady, model, YY)
            #fdxdy = posterior!(paradxdy, model, YY)

            #gradx[k] = -(fx - fdx) / (dx[k]*dxscale[seli])
            #grady[k] = (fx - fdy) / (dx[k]*dxscale[seli])
            #gradxy[k] = -(fx - fdxdy) / sqrt((dx[k]*dxscale[seli])^2 + (dx[k]*dxscale[seli])^2)
            hessdiag[k] = -(2fx - fdx - fdy) / (dx[k]*dxscale[seli])^2
            #hessdiag[k] = -(fx - fdx - fdy + fdxdy) / (dx[k]*dx[k]*dxscale[seli]*dxscale[seli])
        end

        if noisy
            println("Values: $(-hessdiag)")
        end

        hessian[seli, seli] = -(hessdiag[3]+hessdiag[4])/2
        if hessian[seli, seli] < 0
            error("Negative diagonal in Hessian")
        end

        if noisy
            value = hessian[seli, seli]
            println("Value used: $value\n")
        end
    end

    # Now compute off-diagonal elements
    # Make sure that correlations are between -1 and 1
    # errorij contains the index of elements that are invalid
    errorij = Dict{(Int64, Int64), Float64}()

    for i = 1:(nfree-1)
        seli = fpara_free[i]
        for j = (i+1):nfree
            selj = fpara_free[j]

            if noisy
                println("Hessian element: ($seli, $selj)")
            end

            for k = 1:ndx
                paradx = copy(x)
                parady = copy(x)
                paradx[seli] = paradx[seli] + dx[k]*dxscale[seli]
                parady[selj] = parady[selj] - dx[k]*dxscale[selj]
                paradxdy = copy(paradx)
                paradxdy[selj] = paradxdy[selj] - dx[k]*dxscale[selj]

                fx = posterior!(x, model, YY)
                fdx = posterior!(paradx, model, YY)
                fdy = posterior!(parady, model, YY)
                fdxdy = posterior!(paradxdy, model, YY)

                #gradx[k] = -(fx - fdx) / (dx[k]*dxscale[seli])
                #grady[k] = (fx - fdy) / (dx[k]*dxscale[selj])
                #gradxy[k] = -(fx -fdxdy) / sqrt((dx[k]*dxscale[selj])^2 + (dx[k]*dxscale[seli])^2)
                #hessdiag[k] = -(2fx - fdx - fdy) / (dx[k]*dxscale[seli])^2
                hessdiag[k] = -(fx - fdx - fdy + fdxdy) / (dx[k]*dx[k]*dxscale[seli]*dxscale[selj])
            end

            if noisy
                println("Values: $(-hessdiag)")
            end

            hessian[seli, selj] = -(hessdiag[3]+hessdiag[4])/2

            if hessian[seli, selj] == 0 || hessian[selj, selj] == 0
                corrij = 0
            else
                corrij = hessian[seli, selj] / sqrt(hessian[seli, seli]*hessian[selj, selj])
            end

            if corrij < -1 || corrij > 1
                hessian[seli, selj] = 0
                errorij[(seli, selj)] = corrij
            end

            hessian[selj, seli] = hessian[seli, selj]

            if noisy
                value = hessian[seli, selj]
                println("Value used: $value")
                println("Correlation: $corrij")
                println("Number of errors: $(length(errorij))\n")
            end
        end
    end

    stoph = false
    if !isempty(errorij)
        println("Errors: $errorij")
        stoph = true
    end

    return hessian, stoph
end
