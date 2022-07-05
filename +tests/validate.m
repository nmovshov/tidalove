function T = validate(N)
% VALIDATE Test k2 calculation against known cases.

V = [];
names = {};

%% Uniform density
zvec = linspace(1/N, 1, N);
dvec = ones(size(zvec));
k2_calc = lovek2(zvec,dvec);
k2_expect = 1.5;
k2_err = abs(k2_calc - k2_expect)/abs(k2_expect + eps);
row = [k2_expect, k2_calc, k2_err];
names{end+1} = 'Uniform density';
V = [V; row];

%% n=1 polytrope (approximating Jupiter-like planet)
M = 317.8*5.9722e+24;
R = 71492e3;
G = 6.67430e-11;
K = 2*G/pi*R^2;
a = sqrt(2*pi*G/K);
R = pi/a; % should reproduce R
rho_av = 3*M/(4*pi*R^3);
rho_c = (pi^2/3)*rho_av;
r = linspace(1/N,1,N)*R;
rho_exact = rho_c*sin(a*r)./(a*r);
rho_exact(end) = 0; % could end up on the negative side of double prec zero
zvec = flipud(r);
dvec = flipud(rho_exact);
k2_calc = lovek2(zvec,dvec);
k2_expect = (15/pi^2) - 1;
k2_err = abs(k2_calc - k2_expect)/abs(k2_expect + eps);
row = [k2_expect, k2_calc, k2_err];
names{end+1} = 'n=1 poly';
V = [V; row];

%% Tabelize
T = table();
T.expected = V(:,1);
T.calculated = V(:,2);
T.rel_err = V(:,3);
T.Properties.RowNames = names;
end
