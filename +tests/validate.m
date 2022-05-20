function T = validate(N)
% VALIDATE Test k2 calculation against known cases.

M = [];
zvec = linspace(1/N, 1, N);

% Uniform density
dvec = ones(size(zvec));
k2_calc = lovek2(zvec,dvec);
k2_expect = 1.5;
k2_err = abs(k2_calc - k2_expect)/abs(k2_expect + eps);
row = [k2_expect, k2_calc, k2_err];
M = [M; row];


% Tabelize
T = table();
T.expected = M(:,1);
T.calculated = M(:,2);
T.rel_err = M(:,3);
end
