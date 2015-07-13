%% This file generates the variables in gensys_args.mat
% G0, G1, C, PSI, PIE, div
% Beginning with Main.m, it evaluates all code up to but not including the
% call to gensys in dsgesolv.m

cd P:\LSE_2015-FRBNYDSGE_code\cleanCode990

%% Main.m

clear
spec_990
set_paths



%% initializePrograms.m

% Set Random Number Generators
rand('state',12345)
randn('state',12345)

% Initialize Model Specifications
spec

% Load data
[YYall,XXall,ti,nobsall,dlpopall,dlMA_pop,MA_pop,population] = loaddata(nvar,nlags,nant,antlags,psize,zerobound,peachflag,mspec);

% Augment number of observerved variables by number of anticipated policy shocks.
nvar = nvar+nant;

% Define in-sample data
I = find(ti == dates);
YY = YYall(I-stime+1:I,:);
XX = XXall(I-stime+1:I,:);
tiI = ti(I-stime+1:I);
dlpop = dlpopall(I-stime+1:I);
pop_smth = population(I-stime+2+nlags:I+1+nlags);
nobs = size(YY,1);

% Define pre-sample data
YY_p = YYall(1:I-stime,:);
XX_p = XXall(1:I-stime,:);
dlpop_p = dlpopall(1:I-stime);
nobs_p = size(YY_p,1);

% Initialize date-relevant variables
Idate = find(tiI == dates);
tiall    = [tiI(1:end-1);(tiI(end):0.25:(tiI(end)+.25*(qahead+stime-Startdate)))'];
datesall = strcat(num2str(floor(tiall)),'-',num2str(round(1+4*(tiall-floor(tiall)))));

% Transformations associated with level variables
if mspec==557
    YY(1,[1:6 8:end]) = NaN;
    YY(2:end,7)       = NaN;
end

% Load information for prior distribution
prior = priors990();

pmean  = prior(1:npara,1);
pstdd  = prior(1:npara,2);
pshape = prior(1:npara,3);

pshape = pshape.*(1-para_mask);
pmean = pmean.*(1-para_mask)+para_fix.*para_mask;
pstdd = pstdd.*(1-para_mask);

nonpolipar = [1:1:npara];
nonpolipar(polipar) = [];

% Load transformation scheme for parameters
trspec = transp990();
trspec(:,1) = trspec(:,1).*(1-para_mask);



%% dsgesolv.m

valid = 1;
TTT = 0;
RRR = 0;
CCC = 0;

% Assign a value to each parameter
getPara_script

% Pre-allocate matrices


if exist('fzflag','var')
    if fzflag < 1
        disp('fzero failed to converge');
        TTT = [];
        RRR = [];
        valid = -1;
        return
    end
end

% Assign a number to each state

eval(['states',num2str(mspec)]);

% Assign numbers to each equations

eval(['eqs',num2str(mspec)]);

% Express the equilibrium conditions in canonical form

% Total number of equations, states, shocks, and endogenous variables

nstate = n_end+n_exo+n_exp;
neq  = nstate;

G0 = zeros(neq,neq);
G1 = zeros(neq,neq);
C =  zeros(neq,1);
PSI = zeros(neq,nex);
PIE =  zeros(neq,nend);


eval(['eqcond',num2str(mspec)]);

if any(isnan(G0(:))) || any(isnan(G1(:)))
    disp('NaN values in G0 or G1');
    %keyboard;
    valid = -2;
    return
end



%% Write arguments

div = 1 + 1e-6;

save('P:\test_Gensys\gensys_args.mat', 'G0', 'G1', 'C', 'PSI', 'PIE', 'div');