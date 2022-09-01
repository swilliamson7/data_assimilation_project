% script computes the same thing as rossby_wave, but added in the ability
% to plot the actual contours we get from the solution (so a sum of the
% normal modes). we can also include the time-independent Stommel solution

% this is just a test of the timestepping of the rossby wave equation. I'm
% using the analytical solution given in Carl's paper that includes **no
% forcing**, so all results depend on the initial condition. We assume a
% square domain with side length L, and psi is zero along all boundaries. 

% The analytical solution is given by 
%
%  psi(x, y, t) = sum_m sum_n exp(-i sigma_{nm} t) exp(- i beta x/
%                  sigma_{nm}) sin(n pi x) sin(m pi y)
%
% where
%           beta - parameter from the beta-plane approximation for the
%           coriolis force
%           sigma_{nm} - dispersion relation

% we non-dimensionalize the Rossby wave equation so that beta is one and L
% is one.

% We then have the equation for sigma_{nm}:
%        
%  sigma_{nm} = beta / ( 2 (m^2 pi^2 + n^2 pi^2)^(1/2) )
%

dt = 29;
L = 1; 
beta = 1;
f0 = 1e-4;
T = 10;

% deciding which modes to keep
n = 3:5;
m = 4:9;

dx = .01; dy = .01;

% not sure if we actually need to keep this function
f = @(y) f0 + beta * y;

x = 0:dx:L; y = 0:dy:L;
nx = length(x); ny = length(y);
[X, Y] = meshgrid(x,y);

psi = zeros(nx, ny);

N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  

sigma_nm = beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );

% with no dissipation
A = diag( exp(-1i * sigma_nm(:) * dt ) );

state = ones(M*N,1);
mode_coeff = zeros(M*N,T);

b = 1/553;
for j = 1:T
    
    % uncomment to add numerical dissipation
    A = exp(-b * j * 1/0.29) .* A;
    
    state = A * state;
    mode_coeff(:,j) = state;
    
end

% now we use the computed coefficients to find the contour plot of the
% rossby wave 

vec_sigma_nm = sigma_nm(:);
vec_n = Nn(:);
vec_m = Mm(:);

for j = 1:M*N
        
        % no Stommel solution
        
        psi = psi + mode_coeff(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
            .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y);
        
        % with the added stommel solution
        
%         psi = psi + mode_coeff(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
%             .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y) ...
%             + exp(-X) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
        

end

contourf(x, y, real(psi))
colorbar
