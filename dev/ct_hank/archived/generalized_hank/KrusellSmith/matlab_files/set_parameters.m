%% Set Parameters
% Preferences
global ggamma rrho
ggamma = 2;        % coefficient of relative risk aversion
rrho = 0.01;       % rate of time preference

% Production function
global ddelta aalpha
ddelta = .025;     % capital depreciation
aalpha = 1 / 3;    % capital share

% Aggregate shock
global ssigmaTFP rrhoTFP
ssigmaTFP = .007;  % standard deviation of TFP shock
rrhoTFP = .95;     % quarterly autocorrelation of TFP shock

% Idiosyncratic shocks
global z
zz1 = 0;           % unemployed
zz2 = 1;           % employed
z = [zz1,zz2];

% Transition probabilities
global lla
llambda1 = 1 / 2;  % expected duration of unemployment is 2 quarters
llambda2 = (llambda1 / (zz2 * .93 - zz1))*(zz2 - zz2 * .93); % unemployment rate 7%
lla = [llambda1,llambda2];

% Tax system
global mmu ttau
mmu = .15;        % UI replacement rate 15%
ttau = (mmu / zz2) * (lla(2) / lla(1));	     % labor income tax

% Labor supply
global zAvg
zAvg = (lla(1) * z(2) + lla(2) * z(1)) / (lla(1) + lla(2));

%% Approximation Parameters
% Wealth grid
global I amin amax a da aa aaa
I = 100;           % can't be too fine for some reason
amin = 0;
amax = 100;
a = linspace(amin,amax,I)';
da = (amax - amin) / (I - 1);
aa = [a,a];
aaa = reshape(aa,2*I,1);

% Labor productivity
global zz zzz
zz = ones(I,1) * z;
zzz = reshape(zz,2*I,1);

% Idiosyncratic shocks for income
global Aswitch
Aswitch = [-speye(I) * lla(1),speye(I) * lla(1);speye(I) * lla(2),-speye(I) * lla(2)];

% Steady state computations
global rmin rmax r0 maxit crit Delta Ir crit_S
rmin = .0001;       % lower bound for steady state interest rate
rmax = rrho;        % upper bound for steady state interest rate
r0 = .005;          % initial guess for steady state interest rate
maxit = 100;        % maximum iterations on steady state HJB
crit = 1e-6;        % error criterion for steady state value function convergence
Delta = 1e4;        % update size for implicit scheme on steady state HJB
Ir = 100;           % maximum iterations on steady state interest rate
crit_S = 1e-5;      % error criterion for steady state interest rate

% Number of variables in the system
global nVars nEErrors n_v n_g n_p
n_v = 2 * I;
n_g = 2 * I-1 + 1;
n_p = 6;
nVars = n_v + n_g + n_p;
nEErrors = 2 * I;

% Spline parameters
global n_splined
n_knots = 12;
c_power = 7;
x = a';
n_post = 2;	        % This is from two income states
n_prior = 1;
n_splined = n_prior*n_knots*n_post;