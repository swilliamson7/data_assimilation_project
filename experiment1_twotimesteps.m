clc, clear

% This script will run the first experiment: near perfect observations at
% two time steps. 

x0 = [1;0;0;0;0;0];
k = 30;
r = 0.5;
dt = 0.001;
M = 10000;

% known forcing, only applied to position of mass one
B = x0;
q = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

% unknown forcing, again only applied to position of mass one
Gamma = x0;
load('random_forcing_0_1_sd.mat')

% building continuous operators
[Ac, ~, ~] = build_matrices(k, r);

% discrete forward operator 
A = eye(6) + dt .* Ac; 

% create the "data" for the KF
load('noise_0_1_sd.mat')
[A, Ac, Kc, Rc, all_states, eps, k_forcing, ...
    ~, ~] = forward_func(x0, k, r, dt, M, B, q, Gamma, u);

% initial values for the Kalman filter
data = all_states + noise;

x0_KF = x0;

P0 = zeros(6);

E = eye(6);

Q = var(u); 

R = var(noise) .* eye(6);

% known forcing seen by KF
q_KF = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

% create the prediction model, one with no random forcing but correct
% initial conditions, A, B matrices
[~, ~, ~, ~, all_states_pred, eps_pred, k_forcing_pred, kinetic_pred, ...
                potential_pred] = forward_func(x0, k, r, dt, M, B, q_KF, 0.0*Gamma, 0.0.*u);

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

% place to store covariance matrix P
uncertainty = cell(M,1);
uncertainty{1} = P0;

state_old = x0;
P_old = P0;
temp = dt;
for j = 2:M
    
    k_forcing(j) = q_KF(temp);
%     forcing = k_forcing(j) + u(j);
    
    % step forward
    temp_state = A * state_old + B * k_forcing(j);
    temp_P = A * P_old * A' +  (Gamma * Q * Gamma');
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if j == 5000 || j == 7300
        
        state_now = temp_state + K * (data(:, j) - E * temp_state );
        P_now = temp_P - K * E * temp_P;
    else
        state_now = temp_state;
        P_now = temp_P;
        
    end
    
    % calculate and store system energy
    kinetic = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps_KF(j) = kinetic + potential;

    % store new state
    all_states_KF(:,j) = state_now;
    
    % store new uncertainty 
    uncertainty{j} = P_now;

    % cycle states
    state_old = state_now;
    P_old = P_now;
    
    temp = temp + dt;
    
end 

t = 0:M-1;

tiledlayout(3,1);

energy_comparison = nexttile;
plot(energy_comparison, t, eps, t, eps_pred, '--', 'linewidth', 1.5)
ylabel('Energy')
legend('True energy', 'Prediction energy')

kalman_mass1_vel = nexttile;
std_dev = zeros(1,M);
for j = 1:M
    
    std_dev(j) = sqrt(uncertainty{j}(4,4));
    
end
std_dev_plot = errorbar(kalman_mass1_vel, t, all_states_KF(4,:), std_dev);
%alpha = 0.65;   
% Set transparency (undocumented)
%set([std_dev_plot.Bar, std_dev_plot.Line], 'ColorType', 'truecoloralpha', 'ColorData', [std_dev_plot.Line.ColorData(1:3); 255*alpha])
xline([5000, 7300], ':')
ylabel('Displacement')
legend('x_4(t)')


kalman_energy = nexttile;
plot(kalman_energy, t, eps, t, eps_KF, '--', 'linewidth', 1.5)
xline([5000, 7300], ':')
ylabel('Energy')
legend('True energy', 'KF energy')

xlabel('Time step')
