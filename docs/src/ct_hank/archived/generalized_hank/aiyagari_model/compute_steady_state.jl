## Solves the (non-stochastic) steady state of one-asset HANK model

function compute_steady_state(grids::Dict{Symbol, Array{Float64}},
                              params::Dict{Symbol, Float64},
                              approx_params::Dict{Symbol, Float64};
                              Ir::Int64 = 100, maxit_HJB::Int64 = 500,
                              tol_HJB::Float64 = 1e-8, Delta_HJB::Float64 = 1e6,
                              maxit_KFE::Int64 = 1000, tol_KFE::Float64 = 1e-12,
                              Delta_KFE::Float64 = 1e6, niter_hours::Int64 = 10,
                              IterateR::Bool = false, IterateRho::Bool = true,
                              DisplayLev::Int64 = 2)


#Finite difference approximation of the partial derivatives
V = zeros(I,J);
Vaf = zeros(I,J);
Vab = zeros(I,J);
Vzf = zeros(I,J);
Vzb = zeros(I,J);
Vzz = zeros(I,J);
Va_Upwind = zeros(I,J);
c = zeros(I,J);
cf = zeros(I,J);
cb = zeros(I,J);
sf = zeros(I,J);
sb = zeros(I,J);
c0 = zeros(I,J);
If = zeros(I,J);
Ib = zeros(I,J);
I0 = zeros(I,J);
X = zeros(I,J);
Y = zeros(I,J);
Z = zeros(I,J);
A = spzeros(I*J,I*J);
AT = spzeros(I*J,I*J);
b = zeros(I*J);

#CONSTRUCT MATRIX Aswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;          # -σ²/(Δz²) - μ/(Δz)
chi =  s2/(2*dz2);              # σ²/(2*Δz²)
zeta = mu/dz + s2/(2*dz2);      # μ/(Δz) + σ²/(2*Δz²)

updiag = zeros(I*J);
centdiag = zeros(I*J);
lowdiag = zeros(I*J);

for i = 1:I
    centdiag[i]=chi[1]+yy[1];
    for j = 1:J-2
        centdiag[I*j+i] = yy[j+1];
        lowdiag[I*(j-1)+i] = chi[j+1];
        updiag[I*j+i] = zeta[j];
    end
    centdiag[(J-1)*I+i] = yy[J]+zeta[J];
    updiag[(J-1)*I+i] = zeta[J-1];
    lowdiag[I*(J-2)+i] = chi[J];
end

#Add up the upper, center, and lower diagonal into a sparse matrix
Aswitch = spdiagm(centdiag,0,I*J,I*J)+spdiagm(lowdiag[1:I*J-I],-I,I*J,I*J)+spdiagm(updiag[I+1:I*J],I,I*J,I*J);
        # iterate HJB
        for ihjb = 1:maxit_HJB

            # Grab guess
            V = v

            # Create differences matrices. User must define these matrices.
            Vaf, Vab, cf, hf, cb, hb = construct_initial_diff_matrices(V, Vaf, Vab, income,
                                                                       labor, h, h0, zz,
                                                                       profshare,
                                                                       lumptransfer, amax, amin,
                                                                       coefrra, r, dazf,
                                                                       dazb, maxhours)

#Clear updiag and lowdiag
for j = 1:J
   updiag[I*(j-1)+1] = 0;
   lowdiag[I*j] = 0;
end


#----------------------------------------------------
#INITIAL GUESS
r = alpha     * K^(alpha-1) -delta; #interest rates
w = (1-alpha) * K^(alpha);          #wages
v0 = (w*zz + r.*aa).^(1-ga)/(1-ga)/rho;
v = v0;
dist = zeros(maxit);

#-----------------------------------------------------
#MAIN LOOP
for iter=1:maxitK
    println("Main loop iteration ",iter);

    # HAMILTON-JACOBI-BELLMAN EQUATION %
    for n=1:maxit
        V = v;

        # forward difference
        Vaf[1:I-1,:] = (V[2:I,:]-V[1:I-1,:])/da;
        Vaf[I,:] = (w*z + r.*amax).^(-ga); #will never be used, but impose state constraint a<=amax just in case

        # backward difference
        Vab[2:I,:] = (V[2:I,:]-V[1:I-1,:])/da;
        Vab[1,:] = (w*z + r.*amin).^(-ga);  #state constraint boundary condition

        #I_concave = Vab > Vaf; indicator whether value function is concave (problems arise if this is not the case)

        #consumption and savings with forward difference
        cf = Vaf.^(-1/ga);
        sf = w*zz + r.*aa - cf;

        #consumption and savings with backward difference
        cb = Vab.^(-1/ga);
        sb = w*zz + r.*aa - cb;
        #consumption and derivative of value function at steady state
        c0 = w*zz + r.*aa;
        Va0 = c0.^(-ga);

        # dV_upwind makes a choice of forward or backward differences based on
        # the sign of the drift
        If = sf .> 0; #positive drift --> forward difference
        Ib = sb .< 0; #negative drift --> backward difference

        I0 = (1-If-Ib); #at steady state

        #make sure backward difference is used at amax
        #     Ib[I,:] = 1; If[I,:] = 0;
        #STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS sf > 0:
        #already taken care of automatically

        Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; #important to include third term

        c = Va_Upwind.^(-1/ga);
        u = c.^(1-ga)/(1-ga);

        #CONSTRUCT MATRIX A
        X = - min.(sb,0)/da;
        Y = - max.(sf,0)/da + min.(sb,0)/da;
        Z = max.(sf,0)/da;

    for j=1:J
        for i=1:I
            @inbounds centdiag[I*(j-1)+i] = Y[i,j];
            if i<I
                @inbounds updiag[I*(j-1)+i+1] = Z[i,j];
                @inbounds lowdiag[I*(j-1)+i] = X[i+1,j];
            end
        end
    end

        AA=spdiagm(centdiag,0,I*J,I*J)+spdiagm(updiag[2:I*J],1,I*J,I*J)+spdiagm(lowdiag[1:I*J-1],-1,I*J,I*J);

        A = AA + Aswitch;
        B = (1/Delta + rho)*speye(I*J) - A;

        u_stacked = reshape(u,I*J,1);
        V_stacked = reshape(V,I*J,1);

    BLAS.axpy!(I*J,1/Delta,V_stacked,1,u_stacked,1);
        #b = u_stacked + V_stacked/Delta;

        V_stacked = B\u_stacked[:,1];#B\b[:,1]; #SOLVE SYSTEM OF EQUATIONS

        V = reshape(V_stacked,I,J);

        Vchange = V - v;
        v = V;

        dist[n] = maximum(abs.(Vchange));
        if dist[n]<crit
            println("Value Function Converged, Iteration = ",n);
            break
        end
    end
    # FOKKER-PLANCK EQUATION %
    AT = A';
    for i=1:I*J
        b[i] = 0;
    end

    #need to fix one value, otherwise matrix is singular
    i_fix = 1;
    b[i_fix]=.1;
    for j=1:I*J
    AT[i_fix,j]=0;
    end
    AT[i_fix,i_fix]=1;

    #Solve linear system
    gg = AT\b[:,1];
    g_sum = gg'*ones(I*J,1)*da*dz;
    gg = gg./g_sum;

    g = reshape(gg,I,J);

    # Update aggregate capital
    S = sum(g.*a*da*dz);
    println(S)

    if abs(K-S)<critK
        break
    end

    #update prices
    K = relax*K +(1-relax)*S;           #relaxation algorithm (to ensure convergence)
    r = alpha     * K^(alpha-1) -delta; #interest rates
    w = (1-alpha) * K^(alpha);          #wages

    return vars_SS
end
