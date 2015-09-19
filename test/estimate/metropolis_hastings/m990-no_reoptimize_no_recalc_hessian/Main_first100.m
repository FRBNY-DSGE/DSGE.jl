%% Main_first100.m
%% This script runs the estimation step without reoptimizing/recalculating
%% the hessian, for the first 100 draws of Metropolis-Hastings
% 
% Copyright Federal Reserve Bank of New York.  You may reproduce, use, modify,
% make derivative works of, and distribute and this code in whole or in part
% so long as you keep this notice in the documentation associated with any
% distributed works.   Neither the name of the Federal Reserve Bank of New
% York (FRBNY) nor the names of any of the authors may be used to endorse or
% promote works derived from this code without prior written permission.
% Portions of the code attributed to third parties are subject to applicable
% third party licenses and rights.  By your use of this code you accept this
% license and any applicable third party license.  
% 

% CONDITIONS OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT

% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, EXCEPT TO THE EXTENT
% THAT THESE DISCLAIMERS ARE HELD TO BE LEGALLY INVALID.  FRBNY IS NOT, UNDER
% ANY CIRCUMSTANCES, LIABLE TO YOU FOR DAMAGES OF ANY KIND ARISING OUT OF OR
% IN CONNECTION WITH USE OF OR INABILITY TO USE THE CODE, INCLUDING, BUT NOT
% LIMITED TO DIRECT, INDIRECT, INCIDENTAL, CONSEQUENTIAL, PUNITIVE, SPECIAL OR

% TORT OR OTHER LEGAL OR EQUITABLE THEORY, EVEN IF FRBNY HAS BEEN ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGES OR LOSS AND REGARDLESS OF WHETHER SUCH
% DAMAGES OR LOSS IS FORESEEABLE.

%% Initialization
clear
close all

% set paths
path = '~/dsge/cleanCode990/';
dirs = {[path, '']; [path, 'data/']; [path, 'dsgesolv/']; [path, 'estimation/'];...
    [path,'forecast/']; [path, 'initialization/']; [path, 'kalman/'];...
    [path,'toolbox/']; [path, 'plotting/']};
addpath(dirs{:});

%spath = '/home/rceexm08/.julia/v0.3/DSGE/test/estimate/metropolis_hastings/m990-no_reoptimize_no_recalc_hessian/';
%fpath = '/home/rceexm08/.julia/v0.3/DSGE/test/estimate/metropolis_hastings/m990-no_reoptimize_no_recalc_hessian/matlab_tables/';
% sets important variables and flags
spec_990  

% load in pre-computed random vectors/values
randfile = '/data/dsge_data_dir/dsgejl/estimate/save/input_data/rand_save_big.h5'

%randfile ='/home/rceexm08/.julia/v0.3/DSGE/save/m990/input_data/rand_save.h5'
randvecs = hdf5read(randfile,'randvecs');
randvals = hdf5read(randfile,'randvals');

% Make sure we dont recalc hessian or mode
reoptimize = 0
CH = 0
testing=1

% Settings for super small metropolis-hastings
nblocks = 1;
nsim    = 100;
ntimes  = 1;
nburn   = 0;

%spath = ''
if ~reoptimize & ~CH
  %spath = [pwd, '/no_reoptimize_no_recalc_hessian' ];
  spath = [pwd, '/']
elseif ~reoptimize & CH 

end
fpath = [spath, '/matlab_tables'];
  
%% Estimation
% In this stage we draw from the distribution for the parameters. The modal
% parameters as well as the draws of parameters, are outputted in the /save
% folder.
gibb;

%% Write draws to HDF5 file
num1 = npara*nsim*4     % number of bytes per file
fn = [spath, '/params']
fid = fopen(fn, 'r')


theta = []
while(ftell(fid) < num1)
  
  thetadd = fread(fid,[npara,nsim],'single')';
  theta = [theta;thetadd];
  clear thetadd;
  
end;

fclose(fid)

dset_details.Location = spath;
dset_details.Name = 'theta';

hdf5write('metropolis_hastings.h5', dset_details, 'theta', 'WriteMode','overwrite')



%exit;

