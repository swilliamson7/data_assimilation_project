% This script will solve the Rossby wave equation, now with forcing. In
% order to account for forcing we enforce different boundary conditions:
% our system will now be periodic in x and zero on the y boundaries. 

% The analytical solution is still given by 
%
%  psi(x, y, t) = sum_m sum_n exp(-i sigma_{nm} t) exp(- i beta x/
%                  sigma_{nm}) sin(n pi x) sin(m pi y)
%
% where
%           beta - parameter from the beta-plane approximation for the
%           coriolis force
%           sigma_{nm} - dispersion relation

% for the beta-plane approximation ****add to this****

% We then have the equation for sigma_{nm}:
%        
%  sigma_{nm} = beta L / ( 2 (m^2 pi^2 + n^2 pi^2)^(1/2) )
%
% assuming we have not non-dimensionalized and so L =/= 1.

dt = 120;
L = 1; 
beta = 1;
f0 = 1e-4;
T = 10;

dx = .01; dy = .01;

f = @(y) f0 + beta * y;

x = 0:dx:L; y = 0:dy:L;
[X, Y] = meshgrid(x,y);

% deciding which modes to keep
n = 3:5;
m = 4:9;

N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  

sigma_nm = beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );

sigma = @(n, m) beta * L / ( 2 * pi * sqrt(n^2 + m^2));

A = diag( exp(-1i * sigma_nm(:) * dt ) );

% state_nm = @(t, n, m) exp(-1i * sigma(n, m) * t) * ...
%                     exp(-1i * beta * x / sigma(n, m));
           
state = ones(M*N,1);
mode_coeff = zeros(M*N,T);

for j = 1:T
   
    state = A * state;
    mode_coeff(:,j) = state;
    
end

% creating my figure 10 

figure()
tiledlayout(3,1);

% non-dimensional periods

sigma_inverse = nexttile;
plot(1./sigma_nm(:), '*');

sigma_plot = nexttile;
plot(sigma_nm(:), '*');







