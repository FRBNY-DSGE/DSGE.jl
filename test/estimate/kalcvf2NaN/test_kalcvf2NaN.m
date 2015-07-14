%% Main.m

clear
close all
spec_990  % sets important variables and flags
set_paths % adds necessary paths



%% gibb.m

% Step 1: Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Basic Initialization (Reads in data and model specifications)
initializePrograms;
marglh = 'marglh';

% Determine whether 2part estimation or not
if is2part(mspec)
    % If 2-part estimation, we want to pass nant and antlags to likelihood and
    % objective functions. Also want to pass nant to dsgesolv
    args_nant_antlags = {nant, antlags};
    arg_nant = {nant};
else
    % If not 2-part estimation, pass nant and antlags as 0 to likelihood and
    % objective functions. Don't want to pass nant to dsgesolv
    args_nant_antlags = {0, 0};
    arg_nant = {};
end


%  Configure the Metropolis Algorithm
% Setting the jump size for sampling
cc0 = 0.01;
cc = .09;
date_q = (1:1:qahead)';


% Evaluate posterior at a parameter vector, para
MIN = 0;



%% objfcndsge.m
if MIN
    para = trans(para,trspec);
end
para = para.*(1-para_mask)+para_fix.*para_mask;



%% dsgelh.m

% Set up the general classes of 2-part estimation models where
% operations/indexing is similar within a particular class.
class2part;

% Create structure to hold all matrices -- data matrices, state equation matrices, measurement matrices, etc.
[mt, pd] = dsgelh_partition(YY0, YY, nvar, nant, antlags);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: solution to DSGE model - delivers transition equation for the state variables  S_t
%% transition equation: S_t = TC+TTT S_{t-1} +RRR eps_t, where var(eps_t) = QQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the model
[mt(pd).TTT, mt(pd).RRR, mt(pd).CCC, valid] = dsgesolv(mspec,para, mt(end).nant{:});

%% For models using 2part estimation: Get the normal, no ZB model matrices
if any(mspec == class2part_all)

  % Get the starting indices
  [start_ant_state, start_ant_shock, revol_ind] = get_start_ant(mspec, nant);

  [mt(pd-1).TTT, mt(pd-1).RRR, mt(pd-1).CCC]  = ...
    dsgelh_getNoZB(nant, start_ant_state, start_ant_shock, mt(pd).TTT, mt(pd).RRR, mt(pd).CCC,revol_ind{:});
end

if valid < 1;
    pyt = -1E10;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: define the measurement equation: X_t = ZZ*S_t + D + u_t
%% where u_t = eta_t+MM* eps_t with var(eta_t) = EE
%% where var(u_t) = HH = EE+MM QQ MM', cov(eps_t,u_t) = VV = QQ*MM'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some helper functions to create compound matrices by passing the
% period index to access the relevant EE, MM, QQ, RR, VV matrices
makeHH = @(E, M, Q) E + M*Q*M';
makeVV = @(Q, M) Q*M';
makeVVall = @(R, Q, V, H) [[R*Q*R', R*V]; ...
                           [V'*R', H]];


% Get measurement equation matrices set up for all periods that aren't the presample
for p = 2:pd

  % Save measurement equations
  [mt(p).ZZ, mt(p).DD, mt(p).DDcointadd, mt(p).QQ, mt(p).EE, mt(p).MM, retcode] = ...
    feval(['measur',num2str(mspec)], mt(p).TTT, mt(p).RRR,valid,para,mt(p).nvar,nlags,mspec,npara,coint,cointadd,mt(p).nant{:});

  if retcode == 0
      % invalid parameterization
      pyt = -1E10;
      return;
  end;

  mt(p).HH = makeHH(mt(p).EE, mt(p).MM, mt(p).QQ);
  mt(p).VV = makeVV(mt(p).QQ, mt(p).MM);
  mt(p).VVall = makeVVall(mt(p).RRR, mt(p).QQ, mt(p).VV, mt(p).HH);

  if p == 2 && any(mt(p).CCC ~= 0)
    mt(p).DD = mt(p).DD + (mt(p).ZZ)*((eye(size(mt(p).TTT))-mt(p).TTT)\mt(p).CCC);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 3: compute log-likelihood using Kalman filter - written by Iskander
%%         note that Iskander's program assumes a transition equation written as:
%%         S_t = TTT S_{t-1} +eps2_t, where eps2_t = RRReps_t
%%         therefore redefine QQ2 = var(eps2_t) = RRR*QQ*RRR'
%%         and  VV2 = cov(eps2_t,u_u) = RRR*VV
%%         define VVall as the joint variance of the two shocks VVall = var([eps2_t;u_t])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PRESAMPLE

% Solve lyapunov with normal period state matrices (i.e. period 2 matrices)
[A0,P0] = lyap_nonstationary(mspec,para,mt(2).TTT,mt(2).RRR,mt(2).QQ);

% If there is a presample
if ~isempty(mt(1).YY)

  % If number of presample series differ from number of measurement equation
  % series, chop off extra measurement equations and recompute HH, VV, VVall
  % for presample. Otherwise, HH, VV, VVall same in presample as in normal pd
  if size(mt(1).YY,2) ~= size(mt(2).ZZ,1)
    nvar0 = size(mt(1).YY,2);

    % Chop off parts from measurement matrices
    meas_mats = {'ZZ', 'DD', 'EE', 'MM'};
    for mm = 1:length(meas_mats)
      mt(1).(meas_mats{mm}) = mt(2).(meas_mats{mm})(1:nvar0,:);
    end
    mt(1).EE = mt(1).EE(:,1:nvar0);

    % Recompute
    mt(1).HH = makeHH(mt(1).EE, mt(1).MM, mt(2).QQ);
    mt(1).VV = makeVV(mt(2).QQ, mt(1).MM);
    mt(1).VVall = makeVVall(mt(2).RRR, mt(2).QQ, mt(1).VV, mt(1).HH);

  else
    mt(1).ZZ = mt(2).ZZ;
    mt(1).DD = mt(2).DD;
    mt(1).VVall = mt(2).VVall;
  end

  % run Kalman filter on initial observations
  data = (mt(1).YY)';
  lead = 1;
  a = zeros(size(mt(2).TTT, 2), 1);
  F = mt(2).TTT;
  b = mt(1).DD;
  H = mt(1).ZZ;
  var = mt(1).VVall;
  z0 = A0;
  vz0 = P0;
  save('P:\test\kalcvf2NaN_args.mat', 'data', 'lead', 'a', 'F', 'b', 'H', 'var', 'z0', 'vz0')
  
  
  [L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt] = kalcvf2NaN(data, lead, a, F, b, H, var, z0, vz0);
  save('P:\test\kalcvf2NaN_out9.mat', 'L', 'zend', 'Pend', 'pred', 'vpred', 'yprederror', 'ystdprederror', 'rmse', 'rmsd', 'filt', 'vfilt')
  
  [L, zend, Pend, pred, vpred, yprederror, ystdprederror, rmse, rmsd, filt, vfilt] = kalcvf2NaN(data, lead, a, F, b, H, var);
  save('P:\test\kalcvf2NaN_out7.mat', 'L', 'zend', 'Pend', 'pred', 'vpred', 'yprederror', 'ystdprederror', 'rmse', 'rmsd', 'filt', 'vfilt')
  
  
else

  mt(1).zend = A0;
  mt(1).Pend = P0;

end