clc, clear 

% this script can be used to mess around with data, types of data, and
% subsequently the covariance matrix. It will be used to answer questions
% regarding how elements of our covariance matrix change 

x0 = [1;0;0;0;0;0];
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

% create the "data" for the KF
load('noise_0_1_sd.mat')
[A, Ac, Kc, Rc, all_states, eps, k_forcing, ...
    ~, ~] = forward_func(x0, k, r, dt, M, B, q, Gamma, u);
data = all_states + noise;

% initial values for the Kalman filter

x0_KF = x0;

P0 = zeros(6);

E = zeros(6);
E(5,5) = 1;

Q = var(u); 

R = var(noise) .* eye(6);

% determining which steps will contain data for the KF
t = 1:20;
dataset_1 = [3000 + t];
%t = 1:24;
%dataset_2 = [7000 + 125*t]; 
dataset_2 = [3000];
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

state_old = x0;
P_old = P0;
temp = dt;
for j = 2:M
    
    k_forcing(j) = q_KF(temp);
    forcing = k_forcing(j) + u(j);
    
    % step forward
    temp_state = A * state_old + B * k_forcing(j);
    temp_P = A * P_old * A' + (Gamma * Q * Gamma');
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == dataset_1) || sum(j == dataset_2)
        
        state_now = temp_state + K * (E * data(:, j) - E * temp_state);
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

before_obs = uncertainty{2999};
after_obs = uncertainty{3010};

figure()
tiledlayout(3,2);

diag_P_positions = nexttile;
y_1 = [before_obs(1,1), before_obs(2,2), before_obs(3,3)];
y_2 = [after_obs(1,1), after_obs(2,2), after_obs(3,3)];

plot(diag_P_positions, 1:3, y_1, 1:3, y_2, '--', 'linewidth', 1.5)
xticks(1:3)
diag_legend = legend('Before observations', 'After observations');
grid on
ylabel('Diagonal elements of P', 'FontSize', 12) 
title('Elements corresponding to positions')

diag_P_velocities = nexttile;

y_1 = [before_obs(4,4), before_obs(5,5), before_obs(6,6)];
y_2 = [ after_obs(4,4), after_obs(5,5), after_obs(6,6)];
    
plot(diag_P_velocities, 4:6, y_1, 4:6, y_2, '--', 'linewidth', 1.5)
xticks(4:6)
legend('Before observations', 'After observations')
grid on
ylabel('Diagonal elements of P', 'FontSize', 12)
title('Elements corresponding to velocities')


rows_P_pos_before = nexttile;

plot(rows_P_pos_before, 1:3, before_obs(1:6,1:3)', 'linewidth', 1.5)
xticks(1:3)
ylabel('Row elems. of P before obs.', 'FontSize', 12)
row_legend = legend('Row 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'Row 6');
row_legend.Layout.Tile='south';
row_legend.Orientation='horizontal';
grid on

rows_P_vel_before = nexttile;
plot(rows_P_vel_before, 4:6, before_obs(1:6,4:6)', 'linewidth', 1.5)
xticks(4:6)
ylabel('Row elems. of P before obs.', 'FontSize', 12)
%legend('Row 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'Row 6')
grid on


rows_P_pos_after = nexttile;

plot(rows_P_pos_after, 1:3, after_obs(1:6,1:3)', 'linewidth', 1.5)
xticks(1:3)
ylabel('Row elems. of P after obs.', 'FontSize', 12)
ylim([-6,10])
%legend('Row 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'Row 6')
grid on

 
xlabel('Element of x', 'FontSize', 12)

rows_P_vel_after = nexttile;
plot(rows_P_vel_after, 4:6, after_obs(1:6,4:6)', 'linewidth', 1.5)
xticks(4:6)
ylabel('Row elems. of P after obs.', 'FontSize', 12)
%legend('Row 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'Row 6')
grid on


xlabel('Element of x', 'FontSize', 12)
