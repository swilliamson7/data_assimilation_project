clc, clear, close all 

% here we assume our observations are of an average of positions of masses 2
% and 3 and an average of the velocities of masses 2 and 3

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
E = [0 1/2 1/2 0 0 0;
     0 0 0 1/2 1/2 0];

load('noise_for_data_1_sd.mat')
[A, Ac, Kc, Rc, all_states, eps, k_forcing, ...
    kinetic, potential] = forward_func(x0, k, r, dt, M, B, q, Gamma, u);

% observations are averages 
data = E * all_states + 0.01.*noise(1:length(all_states));

% initial values for the Kalman filter
x0_KF = x0;

P0 = zeros(6);

Q = 0.1^2; 

R = 0.01^2 .* eye(2);

% determining which steps will contain data for the KF
t = 1:13;
dataset_1 = 3000 + 300*t;
t = 1:24;
dataset_2 = 7000 + 125*t; 

data_points = [dataset_1, dataset_2];

% forcing seen by KF
q_KF = @(t) 0.5*q(t);

% create the prediction model, one with no random forcing but correct
% initial conditions, A, B matrices, and the same forcing as that seen by
% the Kalman filter 
[~, ~, ~, ~, all_states_pred, eps_pred, k_forcing_pred, kinetic_pred, ...
                potential_pred] = forward_func(x0_KF, k, r, dt, M, B, q_KF, 0.0.*Gamma, 0.0.*u);

% storing the states found by the KF
all_states_KF = zeros(6,M);
all_states_KF(:, 1) = x0_KF;

% place to store energy of the system
eps_KF = zeros(1,M);
kinetic_KF = zeros(1,M);
potential_KF = zeros(1,M);
kinetic_KF(1) = 0.5 * ( x0_KF(4:6)' * x0_KF(4:6) );
potential(1) = 0.5 * (- x0_KF(1:3)' * Kc * x0_KF(1:3) );
eps_KF(1) = kinetic_KF(1) + potential_KF(1);

% place to store uncertainty 
uncertainty = cell(M,1);
uncertainty{1} = P0;

state_old = x0_KF;
P_old = P0;
temp = 0;
for j = 2:M
    
    k_forcing(j) = q_KF(temp);
    forcing = k_forcing(j) + u(j);
    
    % step forward
    temp_state = A * state_old + B * k_forcing(j);
    temp_P = A * P_old * A' + 0.1^2 * (Gamma * Gamma');
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == data_points)
        
        state_now = temp_state + K * (data(:, j) - E * temp_state);
        P_now = temp_P - K * E * temp_P;
    else
        state_now = temp_state;
        P_now = temp_P;
        
    end
    
    % calculate and store system energy
    kinetic1 = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential1 = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps_KF(j) = kinetic1 + potential1;
    kinetic_KF(j) = kinetic1;
    potential_KF(j) = potential1;

    % store new state
    all_states_KF(:,j) = state_now;
    
    % store new uncertainty 
    uncertainty{j} = P_now;

    % cycle states
    state_old = state_now;
    P_old = P_now;
    
    temp = temp + dt;
    
end 

x2_uncertainty = zeros(1, M);


for j = 1:M

    x2_uncertainty(j) = uncertainty{j}(2,2)^(1/2);

end


t = 0:M-1;


tiledlayout(4,1);

nexttile;
plot(t, eps, t, eps_pred, '--', t, eps_KF, '-.', 'linewidth', 1.5)
xline(data_points, ':')
ylabel('Energy')
one = legend('$\mathcal{E}(t)$', '$\tilde{\mathcal{E}}(t, -)$', ...
                '$\tilde{\mathcal{E}}(t)$', 'FontSize', 15);
one.Interpreter = "latex";


nexttile;
plot(t, eps - eps_pred, t, eps - eps_KF, '--', 'linewidth', 1.5)
xline(data_points, ':')
ylabel('Energy')
two = legend('$\mathcal{E}(t) - \tilde{\mathcal{E}}(t, -)$', ...
                '$\mathcal{E}(t) - \tilde{\mathcal{E}}(t)$', 'FontSize', 15);
two.Interpreter = "latex";

% nexttile;
% plot(t, all_states(2,:), t, all_states_pred(2,:), '--', t, all_states_KF(2,:), ...
%     '-.', 'linewidth', 1.5)
% ylabel('Displacement')
% xline(data_points, ':')
% three = legend('$x_2(t)$', '$\tilde{x}_2(t, -)$', '$\tilde{x_2(t)}$', ...
%                 'FontSize', 12);
% three.Interpreter="latex";

nexttile;
plot(t, all_states(2,:), 'linewidth', 1.5)
shadedErrorBar(t, all_states_KF(2,:), x2_uncertainty, 'lineprops','-r','patchSaturation',0.33)
hold on 

ylabel('Displacement')
xline(data_points, ':')
% three = legend('$x_2(t)$', '$\tilde{x}_2(t, -)$', '$\tilde{x_2(t)}$', ...
%                 'FontSize', 12);
three = legend('$x_2(t)$', '$\tilde{x_2}(t)$', ...
                'FontSize', 14.5);
three.Interpreter="latex";

nexttile;
plot(t, all_states(2, :) - all_states_pred(2,:), t, all_states(2,:) ...
            - all_states_KF(2,:), '--', 'linewidth', 1.5)
xline(data_points, ':')
yline(0)
four = legend('$x_2(t) - \tilde{x}_2(t,-)$', '$x_2(t) - \tilde{x}_2(t)$', ...
            'FontSize', 14.5);
ylabel('Displacement')
four.Interpreter = "latex";

xlabel('Time step')
