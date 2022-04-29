function [k2, out] = lovek2(varargin)
%LOVEK2 Tidal Love number k2 from density profile.
%   k2 = LOVEK2() returns ....
%
%   [k2, out] = LOVEK2(..., NAME1=VALUE1, NAME2=VALUE2,...)
%   accepts additional parametrs as NAME/VALUE pairs , and also returns an
%   output struct holding diagnostic values and additional derived quantities.
%
% Inputs, positional
% ------------------
% in1 : type1
%     Description.
% in2 : type2
%     Description.
%
% Inputs, NAME/VALUE pairs
% ------------------------
% tol : scalar, positive, (tol=1e-6)
%     Convergence tolerance for something.
% maxiter : scalar, positive, integer, (maxiter=100)
%     Maximum number of some algorithm iterations.
%
% Outputs
% -------
% k2 : scalar, real
%     Love number k2.
% out : struct
%     A structure holding other quantities calculated in the course of running,
%     for diagnosis etc.
%
% Algorithm
% ---------
% Sterne 1939 as described in Buhler 2016.

%% Input parsing
% Zero inputs case, usage only
if nargin == 0
    print_usage()
    return
end
narginchk(0,inf);

%% Return
k2 = 0; % may as well use the latest...
out.iter = 0;

end

%% Helper functions
function print_usage()
    fprintf('Usage:\n\tlovek2(...)\n')
end
