function [k2, out] = lovek2(zvec, dvec, varargin)
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

%% Input handling
% Zero inputs case, usage only
if nargin == 0
    print_usage()
    return
end
narginchk(2,inf);
validateattributes(zvec,{'numeric'},{'real','finite','positive','vector'})
validateattributes(dvec,{'numeric'},{'real','finite','positive','vector'})
zvec = zvec(:);
dvec = dvec(:);
assert(isequal(size(zvec),size(dvec)))

%% Climb up the radius with Buhler (2016) eq. 2
m = dvec(1)*zvec(1)^3; % starting mass
rhom = m/zvec(1)^3; % starting mean density
eta = 0;
for k=1:length(zvec)-1
    deta = (6 - 6*(dvec(k)/rhom)*(eta + 1) + eta - eta^2)/zvec(k);
    eta = eta + deta*(zvec(k+1) - zvec(k));
    m = m + dvec(k+1)*(zvec(k+1)^3 - zvec(k)^3);
    rhom = m/zvec(k+1)^3;
end

%% Return
k2 = (3 - eta)/(2 + eta);
out.iter = 0;

end

%% Helper functions
function print_usage()
    fprintf('Usage:\n\tlovek2(...)\n')
end
