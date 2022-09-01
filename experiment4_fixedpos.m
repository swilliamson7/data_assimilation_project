% This script will run the first experiment: near perfect observations at
% two time steps. 

x0 = [1;0;2;0;0;0];
k = 30;
r = 0.5;
dt = 0.001;
M = 10000;

% known forcing, only applied to position of mass one
B = x0;
q = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

% unknown forcing, a bunch of random numbers with standard deviation 0.1,
% again only applied to position of mass one
Gamma = x0;
load('random_forcing_0_1_sd.mat')

% building continuous operators
[Ac, Kc, Rc] = build_matrices(k, r);

% discrete forward operator 
A = eye(6) + dt .* Ac; 

% loading the noise that will be added to the data points
load('noise_0_01_sd.mat')

% since we're fixing the position, and thus the velocity, of mass 3 we run
% the forward model here for the data points, rather than with the separate
% function 

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

% running the forward model
state_old = x0;
temp = 0;
for j = 2:M
    
    k_forcing(j) = q(temp);
    forcing(j) = k_forcing(j) + u(j);
    
    % step forward
    state_now = A * state_old + B * k_forcing(j) + B * u(j);
    
    % constraining position and velocity
    state_now(3) = 2;
    state_now(6) = 0;

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

data = all_states + noise;

% now we move on to running the Kalman filter 


% initial values for the Kalman filter, assuming correct initial state
x0_KF = x0;

P0 = zeros(6);

E = eye(6);

Q = 0.1^2; 

R = 0.01^2 .* eye(6);

% determining which steps will contain data for the KF
t = 1:20;
dataset_1 = [3000 + 150*t];

% known forcing seen by KF
q_KF = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

% storing the states found by the KF
all_states_KF = zeros(6,M);
all_states_KF(:, 1) = x0_KF;

% place to store energy of the system
eps_KF = zeros(1,M);
kinetic = zeros(1,M);
potential = zeros(1,M);
kinetic(1) = 0.5 * ( x0(4:6)' * x0(4:6) );
potential(1) = 0.5 * (- x0(1:3)' * Kc * x0(1:3) );
eps_KF(1) = kinetic(1) + potential(1);

% place to store uncertainty 
uncertainty = cell(M,1);
uncertainty{1} = P0;

state_old_KF = [1;0;0;0;0;0];
P_old = P0;
temp = 0;
for j = 2:M
    
    k_forcing(j) = q_KF(temp);
    
    % step forward
    temp_state = A * state_old_KF + B * k_forcing(j);
    temp_P = A * P_old * A' + (Gamma * Q * Gamma');
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == dataset_1) 
        
        state_now_KF = temp_state + K * (E * data(:, j) - E * temp_state );
        P_now = temp_P - K * E * temp_P;
    else
        state_now_KF = temp_state;
        P_now = temp_P;
        
    end
    
    % calculate and store system energy
    kinetic = 0.5 * ( state_now_KF(4:6)' * state_now_KF(4:6) );
    potential = 0.5 * ( - state_now_KF(1:3)' * Kc * state_now_KF(1:3) );
    eps_KF(j) = kinetic + potential;

    % store new state
    all_states_KF(:,j) = state_now_KF;
    
    % store new uncertainty 
    uncertainty{j} = P_now;

    % cycle states
    state_old_KF = state_now_KF;
    P_old = P_now;
    
    temp = temp + dt;
    
end 

t = 0:M-1;

tiledlayout(2,1);

displacement_diff = nexttile;
plot(displacement_diff, t, all_states(3,:), t, all_states_KF(3,:), '--', 'linewidth', 1.5)
for j = 1:M
    
    std_dev(j) = sqrt(uncertainty{j}(5,5));
    
end
%errorbar(displacement_diff, t, all_states_KF(5,:), std_dev);
xline(dataset_1, ':')
ylabel('Displacement')
legend('x_3', 'xkf_3')


energy = nexttile;
plot(energy, t, all_states_KF(3,:) - all_states(3,:), 'linewidth', 1.5)
xline(dataset_1, ':')
yline(0)
ylabel('Displacement')
legend('xkf_3 - x_3')



xlabel('Time step')
