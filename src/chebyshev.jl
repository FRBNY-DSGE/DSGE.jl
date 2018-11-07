## Add code
#push!(LOAD_PATH, "/Users/dchilder/Dropbox/DSGE project/matlab/Huggett/Julia/clusterversion")
#to be able to load module

# A bunch of functions for working with Chebyshev polynomials

function dctTypeI(coeffs)
    # Shameless copy/translation of function described below in package Chebfun, ripped from context
    #COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points
    #of the 2nd kind.
    #   V = COEFFS2VALS(C) returns the values of the polynomial V(i,1) = P(x_i) =
    #   C(1,1)*T_{0}(x_i) + ... + C(N,1)*T_{N-1}(x_i), where the x_i are
    #   2nd-kind Chebyshev nodes.
    #
    #   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
    #   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{0}(x_i)
    #   + C(2,j)*T_{1}(x_i) + ... + C(N,j)*T_{N-1}(x_i).
    #
    # See also VALS2COEFFS, CHEBPTS.

    # Copyright 2017 by The University of Oxford and The Chebfun Developers.
    # See http://www.chebfun.org/ for Chebfun information.

    ################################################################################
    # [Developer Note]: This is equivalent to Discrete Cosine Transform of Type I.
    #
    # [Mathematical reference]: Sections 4.7 and 6.3 Mason & Handscomb, "Chebyshev
    # Polynomials". Chapman & Hall/CRC (2003).
    ################################################################################

    # *Note about symmetries* The code below takes steps to
    # ensure that the following symmetries are enforced:
    # even Chebyshev COEFFS exactly zero ==> VALUES are exactly odd
    # odd Chebychev COEFFS exactly zero ==> VALUES are exactly even
    # These corrections are required because the MATLAB FFT does not
    # guarantee that these symmetries are enforced.

    #Convert to floating point type
    coeffs=float(coeffs)

    # Get the length of the input:
    n = size(coeffs, 1);

    # Trivial case (constant or empty):
    if ( n <= 1 )
        fvalues = coeffs;
        return
    end

    # check for symmetry
    isEvencol = maximum(abs.(coeffs[2:2:end,:]),1).== 0;
    isOddcol = maximum(abs.(coeffs[1:2:end,:]),1).== 0;

    # Scale them by 1/2:
    coeffs[2:n-1,:] = coeffs[2:n-1,:]/2;

    # Mirror the coefficients (to fake a DCT using an FFT):
    tmp = [ coeffs ; coeffs[n-1:-1:2,:] ];

    if ( isreal(coeffs) )
        # Real-valued case:
        fvalues = real(fft(tmp,1));
    elseif ( isreal(1i*coeffs) )
        # Imaginary-valued case:
        fvalues = 1i*real(fft(imag(tmp),1));
    else
        # General case:
        fvalues = fft(tmp,1);
    end

    # Flip and truncate:
    fvalues = fvalues[n:-1:1,:];

    # enforce symmetry
    for j=1:size(fvalues,2)
        if (isEvencol[j])
            fvalues[:,j] = (fvalues[:,j]+flipdim(fvalues[:,j],1))/2;
        elseif (isOddcol[j])
            fvalues[:,j] = (fvalues[:,j]-flipdim(fvalues[:,j],1))/2;
        end
    end

    return fvalues

end #dctTypeI

function c2v1(coeffs)
    #c2v1 is a shameless copy of chebfun function described below to get values at type I Chebyshev points

    #COEFFS2VALS   Convert Chebyshev coefficients to values at Chebyshev points
    #of the 1st kind.
    #   V = COEFFS2VALS(C) returns the values of the polynomial V(i,1) = P(x_i) =
    #   C(1,1)*T_{N-1}(x_i) + C(2,1)*T_{N-2}(x_i) + ... + C(N,1), where the x_i are
    #   1st-kind Chebyshev nodes.
    #
    #   If the input C is an (N+1)xM matrix then V = COEFFS2VALS(C) returns the
    #   (N+1)xM matrix of values V such that V(i,j) = P_j(x_i) = C(1,j)*T_{N-1}(x_i)
    #   + C(2,j)*T_{N-2}(x_i) + ... + C(N,j)
    #
    # See also VALS2COEFFS, CHEBPTS.

    # Copyright 2017 by The University of Oxford and The Chebfun Developers.
    # See http://www.chebfun.org/ for Chebfun information.

    #Convert to floating point type
    coeffs=float(coeffs)

    # Get the length of the input:
    n = size(coeffs, 1);
    m = size(coeffs,2);

    # Trivial case (constant):
    if ( n <= 1 )
        return fvalues = coeffs;
    end

    # check for symmetry
    isEvencol = maximum(abs.(coeffs[2:2:end,:]),1).== 0;
    isOddcol = maximum(abs.(coeffs[1:2:end,:]),1).== 0;

    #Precompute weight vector
    w = (exp.(-1im*[0:2*n-1...]*pi/(2*n))/2);
    w[1] = 2*w[1];
    w[n+1] = 0;
    w[n+2:end] = -w[n+2:end];

    # Mirror the values for FFT:
    c_mirror = [coeffs ; ones(1, m) ; coeffs[end:-1:2,:]];

    # Apply the weight vector:
    c_weight = broadcast(*, c_mirror, w);
    fvaltemp = fft(c_weight,1);

    # Truncate and flip the order:
    fvalues = fvaltemp[n:-1:1,:];

    # Post-process:
    if ( isreal(coeffs) )
        # Real-valued case:
        fvalues = real(fvalues);
    elseif ( isreal(1im*coeffs) )
        # Imaginary-valued case:
        fvalues = 1im*imag(fvalues);
    end


    # enforce symmetry
    for j=1:size(fvalues,2)
        if (isEvencol[j])
            fvalues[:,j] = (fvalues[:,j]+flipdim(fvalues[:,j],1))/2;
        elseif (isOddcol[j])
            fvalues[:,j] = (fvalues[:,j]-flipdim(fvalues[:,j],1))/2;
        end
    end

    return fvalues


end #c2v1

function ccquadwts(n,chebkind)
    #CCQUADWTS   Quadrature weights for Chebyshev points of 2nd kind.
    #   CCQUADWTS(N) returns the N weights for Clenshaw-Curtis quadrature on 2nd-kind
    #   Chebyshev points.
    #
    # See also CHEBPTS, BARYWTS.

    # Copyright 2017 by The University of Oxford and The Chebfun Developers.
    # See http://www.chebfun.org/ for Chebfun information.

    ################################################################################
    # DEVELOPER NOTE:
    # We use a variant of Waldvogel's algorithm [1], due to Nick Hale. (See below)
    # We note this is similar to Greg Von Winkel's approach, which can be found on
    # the MathWorks File Exchange.
    #
    # Let $f(x) = \sum_{k=0}^nc_kT_k(x)$, then\vspace*{-3pt} }
    #   I(f) = v.'*c
    # where
    #   v = \int_{-1}^1T_k(x)dx = { 2/(1-k^2) : k even
    #                             { 0         : k odd
    #     = v'*inv(TT)*f(x) where TT_{j,k} = T_k(x_j)
    #     = (inv(TT)'*v)'*f(x)
    # Therefore
    #   I(f) = w.'f(x) => w = inv(TT).'*v;
    # Here inv(TT).' = inv(TT) is an inverse discrete cosine transform of Type I.
    #
    # Furthermore, since odd entries in v are zero, can compute via FFT without
    # doubling up from N to 2N (though we still need to double up from N/2 to N to
    # facilitate the use of ifft).
    #
    # References:
    #   [1] Joerg Waldvogel, "Fast construction of the Fejer and Clenshaw-Curtis
    #       quadrature rules", BIT Numerical Mathematics 46 (2006), pp 195-202.
    #   [2] Greg von Winckel, "Fast Clenshaw-Curtis Quadrature",
    #       http://www.mathworks.com/matlabcentral/fileexchange/6911, (2005)
    ################################################################################

    if ( n == 0 )                      # Special case (no points!)
        w = [];
    elseif ( n == 1 )                  # Special case (single point)
        w = 2;
    elseif (chebkind==2)                                 # General case, Polynomials of 2nd kind
        tvec = 1 .- [2.0:2.0:(n-1)...].^2
        tk=zeros(size(tvec,1)+1,1);
        tk[1,1]=1;
        tk[2:end,1]=tvec
        c = 2 ./ tk';  # Exact integrals of T_k (even)
        cvec=c[floor(Int,n/2):-1:2]';
        ck=zeros(1,size(c,2)+size(cvec,2))
        ck[1,1:size(c,2)] = c;
        ck[1,size(c,2)+1:end] = cvec; # Mirror for DCT via FFT
        intw = real(ifft(ck,2));                   # Interior weights
        zz=intw[1]/2;                      # Boundary weights
        w=zeros(1,n);
        w[1:size(intw,2)]=intw;
        w[1] = zz;
        w[n] = zz;
    elseif (chebkind==1)                                 # General case, Polynomials of 1st kind
        tvec = 1 .- [2.0:2.0:(n-1)...].^2
        tk=zeros(size(tvec,1)+1,1);
        tk[1,1]=1;
        tk[2:end,1]=tvec
        m = 2 ./ tk';  # Exact integrals of T_k (even)
        # Mirror the vector for the use of ifft:
        if mod(n,2)==1 #n is odd
            cvec=-m[floor(Int,(n+1)/2):-1:2]';
            ck=zeros(1,size(m,2)+size(cvec,2))
            ck[1,1:size(m,2)] = m;
            ck[1,size(m,2)+1:end] = cvec; # Mirror for DCT via FFT
        else #n is even
            cvec=-m[floor(Int,n/2):-1:2]';
            ck=zeros(1,size(m,2)+size(cvec,2)+1)
            ck[1,1:size(m,2)] = m;
            ck[1,size(m,2)+2:end] = cvec;
        end
        v = exp.(1im*[0:n-1...]'*pi/n);      # weight (rotation) vector
        c = ck .* v;                      # Apply the weight vector
        w = real(ifft(c,2));             # Call ifft
    else
        println("Error: Chebkind must be 1 or 2")
    end

    return w

end #ccquadwts


function chebpts(n,lo,hi,chebkind)

    #Function to return grid points and associated Clenshaw Curtis quadrature weights for
    #Chebyshev polynomials
    if ( n == 0 )     # Special case (no points)
        x = [];
    elseif n==1
        x = 0;
    else
        if chebkind==2
            #second kind (the ones with points at the boundary)
            m = n - 1;
            # Points for default domain of [-1,1]
            x = copy(transpose(sin.(pi*(-m:2:m)/(2*m))));  # (Use of sine enforces symmetry.)
        elseif chebkind==1
            x = copy(transpose(sin.(pi*((-n+1:2:n-1)/(2*n))))); # (Use of sine enforces symmetry.)
        else
            println("Error: Chebtype Must be 1 or 2")
        end
    end

    xr=((hi .- lo)/2) * x .+ (hi .+ lo)/2; # Rescaled to new grid

    xwts=((hi-lo)/2)*ccquadwts(n,chebkind); #Rescale Clenshaw Curtis Quadrature Weights


    (cgrid,cwts)=[xr,xwts];

end #chebpts

# cheb2leg
# A function to convert Chebyshev coefficients into orthonormal Legendre polynomial coefficients
# Based on cheb2leg.m in Chebfun for Matlab

function cheb2leg(c_cheb,normalize)
    # Input is array of Chebyshev coefficients,
    # set normalize to true to output coefficients on orthonormal legendre basis

    chebsize=size(c_cheb)
    if chebsize[1]==1
        c_leg=c_cheb
    elseif  chebsize[1] > 1
        #     # Use direct approach:
        c_leg = cheb2leg_direct(c_cheb,normalize);
    end

    ## Matlab version gives option to call approximate algorithm when n large
    ## Not included here: maybe useful feature to add
    # else
    #     # Use fast algorithm:
    #     c_leg = cheb2leg_fast( c_cheb, normalize );
    # end

return c_leg

end #cheb2leg



function cheb2leg_direct(c_cheb, normalize)
    #CHEB2LEG_DIRECT   Convert Cheb to Leg coeffs using the 3-term recurrence.
    N = size(c_cheb,1)                     # Number of columns.
    m = size(c_cheb,2)                   # Number of rows.
    N = N - 1;                          # Degree of polynomial.
    x = cos.(.5*pi*transpose((0:2*N))/N);         # 2*N+1 Chebyshev grid (reversed order).
    temparg=zeros(2*N+1,m);
    temparg[1:N+1,1:m]=c_cheb;
    ff = dctTypeI(temparg);                           #Ordering and normalization differ from Matlab version
    f=flipdim(ff,1);                                   # Values on 2*N+1 Chebyshev grid.
    w = copy(transpose(ccquadwts(2*N+1,2)));     # Clenshaw-Curtis quadrature weights.
    Pm2 = 1;
    Pm1 = x;                   # Initialize.
    L = zeros(2*N+1, N+1);              # Vandermonde matrix.
    L[:,1] = 1+0*x;          #P_0
    L[:,2] = x;              #P_1.
    for k = 1:N-1 # Recurrence relation:
        P = (2-1/(k+1))*Pm1 .* x - (1-1/(k+1))*Pm2;
        Pm2 = Pm1;
        Pm1 = P;
        L[:,2+k] = P;
    end
    scale = (2*[0:N...]+1)/2;            # Scaling in coefficients [NIST, (18.3.1)].
    c_leg = broadcast(*, transpose(L)*(broadcast(*, f ,w)), scale); # Legendre coeffs.

    # Normalize:
    if ( normalize )
        c_leg  = broadcast(*, c_leg, 1 ./ sqrt.([0:N...]+.5) );
    end

    return c_leg

end #cheb2leg_direct


function bary(x, fvals, xk)
    #BARY   Barycentric interpolation formula.

    # xk are (2nd kind Chebyshev) grid points, fvals are values at grid points,
    # x are points where you want to interpolate function (using Chebyshev polynomial)

    #   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
    #   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
    #   Note that XK and VK should be column vectors, and FVALS, XK, and VK should
    #   have the same length.
    #
    #   BARY(X, FVALS) assumes XK are the 2nd-kind Chebyshev points and VK are the
    #   corresponding barycentric weights.
    #
    #   If size(FVALS, 2) > 1 then BARY(X, FVALS) returns values in the form
    #   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
    #


    # Copyright 2017 by The University of Oxford and The Chebfun Developers.
    # See http://www.chebfun.org/ for Chebfun information.

    # Parse inputs:
    x=x[:,:]; # Represent as matrix

    n = size(fvals,1);
    m = size(fvals,2);
    sizex = size(x);
    ndimsx = ndims(x);

    if ( (m > 1) && (ndimsx > 2) )
        error("Too many Dimensions")
    end

    # Omit default behavior: always ask for grids and weights
    # # Default to Chebyshev nodes and barycentric weights:
    # if ( nargin < 3 )
    #     xk = chebtech2.chebpts(n);
    # end
#
# if ( nargin < 4 )
#     vk = chebtech2.barywts(n);
# end

#Use Barycentric Weights corresponding to second kind Chebyshev

vk = barywts(n);

# Check that input is a column vector:
if ( (ndimsx > 2) || (sizex[2] > 1) )
    x = x[:];
end

# The function is a constant.
if ( n == 1 )
    return fx = repeat(fvals, size(x,1), 1);
end

# The function is NaN.
if ( any(isnan.(fvals)) )
    return    fx = NaN(size(x,1), m);
end

# The main loop:
if ( prod(size(x)) < 4*size(xk,1) )  # Loop over evaluation points
    # Note: The value "4" here was detemined experimentally.

    # Initialise return value:
    fx = zeros(size(x, 1), m);

    # Loop:
    for j = 1:prod(size(x))
        xvar = vk ./ (x[j]*ones(size(xk)) - xk);
        fx[j,:] = (xvar'*fvals) ./ sum(xvar);
    end
else                            # Loop over barycentric nodes
    # Initialise:
    numerat = zeros(size(x, 1), m);
    denom = numerat;

    # Loop:
    for j = 1:size(xk,1)
        tmp = (vk[j] ./ (x - xk[j]));
        numerat = numerat + broadcast(*, tmp, fvals[j,:]);
        denom = broadcast(+, denom, tmp);
    end
    fx = numerat ./ denom;
end

# Try to clean up NaNs:

nanind=find(isnan,fx[:,1]);

for i = 1:size(nanind,1)
    k=nanind[i]; #Place where there is a NaN
    index = findfirst(xk,x[k]);    # Find the corresponding node
    if  index!=0
        fx[k,:] = fvals[index,:];   # Set entry/entries to the correct value
    end
end

# Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex[2] > 1) ) )
    fx = reshape(fx, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex[2] > 1) ) )
    fx = reshape(fx, sizex[1], m*prod(size(x))/sizex[1]);
end

return fx

end #bary

function barywts(n)
    #returns Barycentric weights for interpolation using Chebyshev polynomials of the second kind

    #BARYWTS   Barycentric weights for Chebyshev points of 2nd kind.
    #   BARYWTS(N) returns the N barycentric weights for polynomial interpolation on
    #   a Chebyshev grid of the 2nd kind. The weights are normalised so that they
    #   have infinity norm equal to 1 and the final entry is positive.
    #
    # See also BARY, CHEBPTS.

    # Copyright 2017 by The University of Oxford and The Chebfun Developers.
    # See http://www.chebfun.org/ for Chebfun information.

    ################################################################################
    # See Thm. 5.2 of Trefethen, Approximation Theory and Approximation Practice,
    # SIAM, 2013 for more information.
    ################################################################################

    if ( n == 0 )                      # Special case (no points)
        v = [];
    elseif ( n == 1 )                  # Special case (single point)
        v = 1;
    else                               # General case
        v = [ones(n-1,1) ; .5];        # Note v(end) is positive.
        v[end-1:-2:1] = -1;
        v[1] = .5*v[1];
    end

    return v

end #barywts

# ChebDerivMat
# Create matrix corresponding to collocation representation of differentiation
# based on recurrence relation for derivatives of Chebyshev polynomials at second kind chebyshev points (Tn(x))
# Returns matrix with proper sclaing, but

function ChebDerivMat(n,lo,hi)


    D = zeros(n,n) #Initialize recurrence relation
    # D is coefficients to coefficients transform
    if n>1
        D[1,2]=1.0 #Polynomial of order 1 has derivative
        if n>2
            D[2,3]=4
            if n>3
                for j=4:n
                    D[j-1,j]=2*(j-1)
                    D[:,j]+=((j-1)/(j-3))*D[:,j-2]
                end
            end
        end
    end

    Minv=dctTypeI(eye(n))
    M=pinv(Minv)

    Deriv=(2/(hi-lo))*Minv*D*M #Grid values to grid values transform

    return Deriv, D

end #ChebDerivMat
