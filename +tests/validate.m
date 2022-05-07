function tf = uniform()
% UNIFORM Test k2 with uniform density.

N = 100;
zvec = linspace(1/N, 1, N);
dvec = ones(size(zvec));
k2 = lovek2(zvec,dvec);

T = table(N,k2);
disp(T)

tf = true;
end
