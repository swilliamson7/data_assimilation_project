% this function will apply the Kalman filter to incorporate data into the
% model. It needs to be passed
%           data - a vector containing the data points
%           x0 - initial condition
%           data_steps1, data_steps2 - where to incorporate data
%           k - spring constant
%           r - dissipation constant 
%           dt - time step size
%           M - number of steps to take 
%           B - distributes known forcing
%           Gamma - distributes unknown forcing
%           q - deterministic forcing function
%           u - random forcing 
%           E - distributes the data points 
%           R - covariance of noise in data 
%           P0 - initial uncertainty
% It will return 
%           all_states_KF - the results of the forward run with data
%           incorporated (if it exists)
%           eps_KF - the energy of the system 
%           uncertainty - the uncertainty of the predicition

function [all_states_KF, eps_KF, uncertainty] = KF_func(x0, k, r, dt, M, B, q, Gamma, u, ...
                                             P0, E, R, data, data_steps1, data_steps2)

[Ac, Kc, Rc] = build_matrices(k,r);

% build discrete forward operator 

A = eye(6) + dt.*Ac;

% place to store state vector
all_states_KF = zeros(6,M);
all_states_KF(:,1) = x0;

% place to store energy of the system
eps_KF = zeros(1,M);
kinetic = zeros(1,M);
potential = zeros(1,M);
kinetic(1) = 0.5 * ( x0(4:6)' * x0(4:6) );
potential(1) = 0.5 * (- x0(1:3)' * Kc * x0(1:3) );
eps_KF(1) = kinetic(1) + potential(1);

% place to store forcing in the system 
forcing = zeros(1,M); 
k_forcing = zeros(1,M);

% place to store uncertainty 
uncertainty = cell(M,1);
uncertainty{1} = P0;

state_old = x0;
P_old = P0;
temp = 0;
for j = 2:M
    
    k_forcing(j) = q(temp);
    forcing(j) = k_forcing(j) + u(j);
    
    % step forward
    temp_state = A * state_old + B * k_forcing(j);
    temp_P = A * P_old * A' + 0.1^2 * (Gamma * Gamma');
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if ismember(j, data_steps1) || ismember(j, data_steps2)
        
        state_now = temp_state + K * (E * data(j) - E * temp_state );
        P_now = temp_P - K * E * temp_P;
    else
        state_now = temp_state;
        P_now = temp_P;
        
    end
    
    % calculate and store system energy
    kinetic(j) = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential(j) = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps_KF(j) = kinetic(j) + potential(j);

    % store new state
    all_states(:,j) = state_now;
    
    % store new uncertainty 
    uncertainty{j} = P_now;

    % cycle states
    state_old = state_now;
    P_old = P_now;
    
    temp = temp + dt;
    
end 

end 
