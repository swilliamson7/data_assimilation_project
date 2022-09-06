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

% building continuous operators
[Ac, Kc, Rc] = build_matrices(k, r);

% discrete forward operator 
A = eye(6) + dt .* Ac; 

% create the "data" for the KF
load('noise_0_01_sd.mat')
[A, Ac, Kc, Rc, all_states, eps, k_forcing, ...
    ~, ~] = forward_func(x0, k, r, dt, M, B, q, Gamma, u);
data = all_states + noise;

% initial values for the Kalman filter
x0_KF = x0;

P0 = zeros(6);

E = [0 1/2 1/2 0 0 0;
     0 0 0 1/2 1/2 0];

Q = 0.1^2; 

R = 0.01^2 .* eye(2);

% determining which steps will contain data for the KF
t = 1:13;
dataset_1 = [3000 + 300*t];
t = 1:24;
dataset_2 = [7000 + 125*t]; 

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
temp = 0;
for j = 2:M
    
    k_forcing(j) = q_KF(temp);
    forcing = k_forcing(j) + u(j);
    
    % step forward
    temp_state = A * state_old + B * k_forcing(j);
    temp_P = A * P_old * A' + 0.1^2 * (Gamma * Gamma');
    
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

tiledlayout(3,1);

displacement_diff = nexttile;
plot(displacement_diff, t, all_states_KF(5,:) - all_states(5,:), 'linewidth', 1.5)
for j = 1:M
    
    std_dev(j) = sqrt(uncertainty{j}(5,5));
    
end
%errorbar(displacement_diff, t, all_states_KF(5,:), std_dev);
xline([dataset_1, dataset_2], ':')
yline(0)
ylabel('Displacement')
legend('xkf_5 - x_5')

energy = nexttile;
plot(energy, t, eps, t, eps_KF, '--', 'linewidth', 1.5)
xline([dataset_1, dataset_2], ':')
ylabel('Energy')
legend('True energy', 'KF energy')

kalman_energy = nexttile;
plot(kalman_energy, t, eps, t, eps_KF, '--', 'linewidth', 1.5)
xline([dataset_1, dataset_2], ':')
ylabel('Energy')
legend('True energy', 'KF energy')
xlim([5000, 8000])

xlabel('Time step')
