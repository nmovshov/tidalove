function k2 = lovek2(zvec, dvec, odeint)
%LOVEK2 Tidal Love number k2 from density profile.
%   k2 = LOVEK2() returns ....
%
% Inputs, positional
% ------------------
% zvec : real finite positive vector
%     Normalized radii starting near (not at) center.
% dvec : real finite nonnegative vector
%     Density at corresponding radius.
% odeint : ['euler', {'trapz'}]
%     Ode integrator type.
%
% Outputs
% -------
% k2 : scalar, real
%     Love number k2.
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
if nargin < 3 || isempty(odeint), odeint = 'trapz'; end
validateattributes(zvec,{'numeric'},{'real','finite','positive','vector'})
validateattributes(dvec,{'numeric'},{'real','finite','nonnegative','vector'})
odeint = validatestring(odeint,{'euler','trapz'});
zvec = zvec(:);
dvec = dvec(:);
assert(isequal(size(zvec),size(dvec)))

%% Climb up the radius with Buhler (2016) eq. 2
m = dvec(1)*zvec(1)^3; % starting mass
rhom = m/zvec(1)^3; % starting mean density
eta = 0;
if odeint == "euler"
    for k=1:length(zvec)-1
        deta = (6 - 6*(dvec(k)/rhom)*(eta + 1) + eta - eta^2)/zvec(k);
        eta = eta + deta*(zvec(k+1) - zvec(k));
        m = m + dvec(k+1)*(zvec(k+1)^3 - zvec(k)^3);
        rhom = m/zvec(k+1)^3;
    end
else
    for k=1:length(zvec) - 1
        s1 = (6 - 6*(dvec(k)/rhom)*(eta + 1) + eta - eta^2)/zvec(k);
        zhalf = zvec(k) + 0.5*(zvec(k+1) - zvec(k));
        dhalf = dvec(k) + 0.5*(dvec(k+1) - dvec(k));
        mhalf = m + dhalf*(zhalf^3 - zvec(k)^3);
        rhalf = mhalf/zhalf^3;
        ehalf = eta + s1*(zhalf - zvec(k));
        s2 = (6 - 6*(dhalf/rhalf)*(ehalf + 1) + ehalf - ehalf^2)/zhalf;
        eta = eta + s2*(zvec(k+1) - zvec(k));
        m = mhalf + dvec(k+1)*(zvec(k+1)^3 - zhalf^3);
        rhom = m/zvec(k+1)^3;
    end
end

%% Return
k2 = (3 - eta)/(2 + eta);
out = [];

end

%% Helper functions
function print_usage()
    fprintf('Usage:\n\tlovek2(zvec, dvec, {odeint=''trapz''})\n')
end
