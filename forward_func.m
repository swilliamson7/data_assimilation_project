% This function will run the forward problem and return the displacements
% at every time step. To run it needs to be passed 
%       dt - time step
%       M - total number of steps to take 
%       x0 - initial state vector
%       k - spring constant
%       r - dissipation constant 
%       B - force distributing matrix
%       q - deterministic forcing function
%       Gamma - matrix that distributes the unknown forcing
%       u - unknown forcing values
% The function returns 
%       A - the full, discrete forward operator 
%       Ac - continuous forward operator
%       Kc, Rc - submatrices of Ac
%       all_states - 6 X M matrix of all states 
%       eps - energy of the system 

function [A, Ac, Kc, Rc, all_states, eps, k_forcing, kinetic, ...
                potential] = forward_func(x0, k, r, dt, M, B, q, Gamma, u) 

[Ac, Kc, Rc] = build_matrices(k,r);

% build discrete forward operator 

A = eye(6) + dt.*Ac;

% place to store state vector
all_states = zeros(6,M);
all_states(:,1) = x0;

% place to store energy of the system
eps = zeros(1,M);
kinetic = zeros(1,M);
potential = zeros(1,M);
kinetic(1) = 0.5 * ( x0(4:6)' * x0(4:6) );
potential(1) = 0.5 * (- x0(1:3)' * Kc * x0(1:3) );
eps(1) = kinetic(1) + potential(1);

% place to store forcing in the system 
forcing = zeros(1,M); 
k_forcing = zeros(1,M);

state_old = x0;
temp = 0;
for j = 2:M
    
    k_forcing(j) = q(temp);
    forcing(j) = k_forcing(j) + u(j);
    
    % step forward
    state_now = A * state_old + B * k_forcing(j) + B * u(j);

    % calculate and store system energy
    kinetic(j) = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential(j) = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps(j) = kinetic(j) + potential(j);

    % store new state
    all_states(:,j) = state_now;

    % cycle states
    state_old = state_now;
    
    temp = temp + dt;
    
end 

end 
