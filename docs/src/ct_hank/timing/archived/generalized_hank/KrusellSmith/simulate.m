function [values,vtime] = simulate(g1,impact,T,N,shocks,method,blowup,subset)
% Given linear dynamics
%        dx = g1*x*dt + impact * dZ
%    simulates the values of x for T time period with N steps with
%    realization of dZ given by shocks
%
% by SeHyoun Ahn, March 2017
%
% REFERENCE: Ahn, SeHyoun, Greg Kaplan, Benjamin Moll, Thomas Winberry, and
%    Christian Wolf. "When Inequality Matters for Macro and Macro Matters
%    for Inequality."
%
% PARAMETERS:
%    dx = g1*x*dt + impact*dZ
%    g1 = (m x m) matrix for dynamics, usually an output of
%         <schur_solver.m>
%    impact = (m x n) matrix of impact of shocks on values, usually an
%             output of <schur_solver.m>
%    T = length of time simulation
%    N = number of steps taken for the simulation
%    shocks = (n x N) matrix of shock realization
%    method = {'implicit','explicit'} updating method for linear system
%    blowup = (optional) translation back to full system if the dynamics is
%             given for a reduced system
%    subset = (optional) subset out variables if only part of the system is
%             needed
%
% OUTPUTS:
%    values = (m x N) matrix of simulated values (subset of the values if
%             optional subset parameter is used)
%    vtime = (N x 1) vector of time values
%
% EXAMPLES:
%    A = diag(-linspace(1,11,10));     % simple stable dynamics
%    impact = linspace(1,11,10)';      % simple impact of shocks
%    T = 5;
%    N = 50;
%    shocks = zeros(1,N); shocks(1) = 1;
%    [vtime,values] = simulate(A,impact,T,N,shocks,'implicit');
%    plot(vtime,values);
%
%    For example with blowup and subset, see Krusell-Smith case example
%       available at < > 
%
% SYNTAX:
% [values,vtime] = simulate(g1,impact,T,N,shocks,method,blowup,subset)

vtime = linspace(0,T,N);
dt = vtime(2)-vtime(1);

% Preallocation
[nvars,~] = size(g1);
values = zeros(nvars,N);

if (method == 'implicit')
    % if N is small, might be faster just to do backslash instead. To use
    %    backslash method, just uncomment line 51/54 and comment line 50/53
    gg1 = inv(speye(size(g1)) - g1*dt);
    % gg1 = speye(size(g1))-g1*dt;
    for n = 1 : N
        values(:,n+1) = gg1*(values(:,n) + (dt^(1/2))*impact*shocks(:,n));
        % values(:,n+1) = gg1\(values(:,n) + (dt^(1/2))*impact*shocks(:,n));
    end
    values = values(:,2:end);
    if nargin > 6
        values = blowup*values;
        if nargin > 7
            values = values(subset,:);
        end
    end
elseif (method == 'explicit')
    gg1 = speye(size(g1))+g1*dt;
    for n = 1 : N
        values(:,n+1) = gg1*(values(:,n) + (dt^(1/2))*impact*shocks(:,n));	
    end
    values = values(:,2:end);
    if nargin > 6
        values = blowup*values;
        if nargin > 7
            values = values(subset,:);
        end
    end
else
    error('<simulate>: Unknown method');
end
