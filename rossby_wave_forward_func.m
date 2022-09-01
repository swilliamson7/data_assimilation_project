% This function will run the forward steps of the Rossby wave system and
% return the coeffecients of the normal mode decomposition. To run it needs 
% to be passed:
%       x0 - initial value, so the initial 18 x 1 vector of coefficients
%       (dependent on the initial condition chosen for psi)
%       dt - time step
%       T - total number of (time) steps to take 
%       n, m - vectors with the values of the modes we'd like to keep
%       beta - beta-plane parameter
%       L - length of the square domain
% The function returns 
%       A - the forward operator
%       sigma_nm - the dissipation relation, given as a length(m) x
%       length(n) matrix, increasing n across rows and increasing m
%       across columns
%       all_states - 18 x M matrix with all the mode coefficients
%       computed
%       energy - energy of the system at every step
%       

function [A, sigma_nm, all_states, energy] = rossby_wave_forward_func(x0, dt, T, n, m, beta, L) 

N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  

sigma_nm = -beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );

A = diag( exp(-1i * sigma_nm(:) * dt ) );
        
state = x0;
all_states = zeros(M*N,T+1);
all_states(:,1) = x0;

energy = zeros(1, T+1);
energy(1) = sum(abs(all_states(:,1)).^2);

for j = 2:T+1
   
    state = A * state;
    all_states(:,j) = state;
    
    energy(j) = sum(abs(state).^2);
    
end

end 