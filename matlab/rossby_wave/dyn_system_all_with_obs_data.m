% This script will compute the coefficients of the modes in our simplified
% dynamical system, not a true solution to the Rossby wave equation. To do 
% this we use the Kalman filter for data assimilation and an RTS smoother
% for smoothing. The difference between this and the script dyn_sys_KF_all.m
% is that here we actually build an operator E capable of moving between wave
% observations and coefficients of our modes, allowing us to 
% examine properties of this matrix later. 
% 
% We define 
%
%      psi(t, x, y) = \sum_n \sum_m exp(-i sigma_nm t - i beta
%           x/sigma_nm) c_nm sin(n pi x) sin(m pi y)
%
% where
%
%      sigma_nm = ( beta/2 ) (1/( pi sqrt(n^2 + m^2) ) )
%
% is our non-dimensional dispersion relation and c_nm are coefficients
% dependent on psi(0, x, y).
%
% We assume the KF has full knowledge of A, the time-stepping matrix, and
% no knowledge of the forcing at every step. We'll provide near perfect
% data at certain timesteps. 
%
% I'm opting to keep this script self-contained, since this will not be a
% true solution of the Rossby wave equation and instead is a simple
% dynamical system that we consider. 

clc, clear, close all

%%

L = 1;
dx = .01; dy = .01;
x = 0:dx:L; y = 0:dy:L;
[X, Y] = meshgrid(x,y);

T = 2000;    % how many steps to take 
dt = 29;     % time step
L = 1;       % length of square basin, 1 b/c non-dimensionalized
b = 1.8e-3;
a = 6.3e6;   % radius of Earth, kilometers
theta = 30;  % latitude we use for computation of beta and f 
Ra = 0.2;    % Stommel parameter, from Carl's code 

% Earth's rotation rate
Omega = 7.2921e-5;         % rads / sec 

% determining beta and f at our given latitude
beta_dim = 2 * Omega * cosd(theta) / a;
f = 2 * Omega * sind(theta);

beta = beta_dim * a / f;

% the choice of n, m determines which of the modes we keep 
n = 3:5;
m = 4:9;
N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  
vec_n = Nn(:);
vec_m = Mm(:);

% initial condition for the system
x0 = (1./(vec_n.^2 + vec_m.^2));

% dispersion relation, given as an N by M matrix, with increases in m
% happening across rows and increases in n happening across columns 
sigma_nm = beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );
vec_sigma_nm = sigma_nm(:);

% creating a function that computess a specific sigma_nm value
sigma_func = @(n,m) -beta * L / (2 * pi * sqrt(n^2 + m^2) );

% building the forward operator for stepping the coefficients 
A = exp(-b) .* diag( exp(-1i * sigma_nm(:) * dt) );    

% random forcing
load('forcing_1_sd.mat');
q = 0.002 .* forcing;

% pre-allocating storage for quantities we want to store

% storage for the "true" states and "energy"
state = x0;
all_states = zeros(M*N,T+1);
all_states(:,1) = x0;

energy = zeros(1, T+1);
energy(1) = sum(abs(all_states(:,1)).^2);

% storage for the KF computed states and energy
x0_KF = 1.2.*x0;
all_states_KF = zeros(M*N, T+1);
all_states_KF(:,1) = x0_KF;        % KF has perfect initial condition

energy_KF = zeros(1, T+1);
energy_KF(1) = sum(abs(all_states_KF(:,1)).^2);

% storage for the predicted states and energy
all_states_pred = zeros(M*N, T+1);
all_states_pred(:,1) = x0_KF;

energy_pred = zeros(1, T+1);
energy_pred(1) = energy_KF(1);

% noise to add to the data points to make them a little fuzzy

load('noise_1_sd.mat');
noise = 0.001.*noise;

%%

% Running the forward model to create the artificial data
for j = 2:T+1
    
    state = A * state + q(:, j);

    all_states(:,j) = state;
    
    energy(j) = sum(abs(state).^2);

end

% now we want to create data from the psi values we just computed. This 
% requires us to build the matrix E differently than in our mass spring 
% oscillator system

x_vals = [1/8:1/8:7/8, 1/8:1/8:7/8];
y_vals = [1/3.*ones(size(1/8:1/8:7/8)), 2/3.*ones(size(1/8:1/8:7/8))];

E = zeros(length(x_vals), M*N);

for j = 1:length(x_vals)
    for k = 1:M*N

        E(j,k) = exp(-1i * beta * x_vals(j) / vec_sigma_nm(k) ) ...
            * sin(vec_n(k) * pi * x_vals(j)) * sin(vec_m(k) * pi * y_vals(j));

    end
end

% % adding noise here to make the data a little fuzzy
data = E * all_states + noise(1:14, 1:T+1); 

%%

% Lastly we want to run the Kalman filter and reconstruct the coefficients
% of the normal mode decomposition 

q_KF = 0.5 .* q;

% covariance matrix
uncertainty = cell(T+1,1);
P0 = cov( (x0 - x0_KF) * (x0 - x0_KF)' );
uncertainty{1} = P0;

uncertainty_before_data = cell(T+1, 1);

% matrix of covariance of noise in the data
R = var(noise(:)) .* eye(length(x_vals));

% covariance of random forcing
Q = var(q_KF(:)) .* eye(N * M);
Gamma = eye(N*M); % we apply forcing to all entries of the state vector

% timesteps where we'll incorporate data
data_steps = [400:50:900, 925:25:1300];
%data_steps = 1:T+1;

% running the Kalman filter and storing the prediction values as we go

state_old = x0_KF;
P_old = P0;
temp = dt;

state_pred = x0_KF;
for j = 2:T+1
    
    state_pred = A * state_pred + q_KF(:, j);
    
    all_states_pred(:, j) = state_pred;
    energy_pred(j) = sum( abs(state_pred).^2 );

end

states_before_data = zeros(M * N, T+1);
states_before_data(:, 1) = x0_KF;

for j = 2:T+1

    temp_state = A * state_old + q_KF(:, j);
    temp_P = A * P_old * A' + Gamma * Q * Gamma'; 
    
    states_before_data(:, j) = temp_state;
    uncertainty_before_data{j} = temp_P;
    
    % calculate Kalman gain matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == data_steps)
        
        state_now = temp_state + K * (data(:, j) - E * temp_state );
        P_now = temp_P - K * E * temp_P;

    else

        state_now = temp_state;
        P_now = temp_P;
        
    end
    
    % computing "energy" 
    energy_KF(j) = sum( abs(state_now).^2 ); 
    
    % store new state
    all_states_KF(:,j) = state_now;
    
    % store new uncertainty 
    uncertainty{j} = P_now;

    % cycle states
    state_old = state_now;
    P_old = P_now;
    
    temp = temp + dt;
    
end 

%%

% Now that the Kalman filter has been used to the final time and all
% available data has been incorporated, we want to run the RTS smoother
% backwards to improve the prediction 

% storage for the smoothed states and energy
all_states_RTS = zeros(M*N, T+1);
energy_RTS = zeros(1, T+1);

uncertainty_RTS = cell(T+1,1);

% the initial value for the smoother is just what was computed with the KF
uncertainty_RTS{end} = uncertainty{end};
all_states_RTS(:, end) = all_states_KF(:, end);

forcing_RTS = zeros(M*N,T+1);

Q_new = Q;

for j = T:-1:1

    L = uncertainty{j} * A' * uncertainty_before_data{j+1}^(-1);
    M_RTS = Q_new * Gamma' * uncertainty_before_data{j+1}^(-1);

    all_states_RTS(:, j) = all_states_KF(:, j) + ...
                L * (all_states_RTS(:, j+1) - states_before_data(:, j+1));
    uncertainty_RTS{j} = uncertainty{j} + L * (uncertainty_RTS{j+1} ...
                - uncertainty_before_data{j+1}) * L';

%     forcing_RTS(:, j) = 0.5 .* q(j) + M * (all_states_RTS(:, j+1) - ...
%                                             states_before_data(:, j+1));
    Q_new = Q_new + M_RTS * (uncertainty_RTS{j+1} - uncertainty_before_data{j+1}) ...
                                * M_RTS';

end

% temp = (A * uncertainty{1} * A' + Gamma * Q * Gamma')^(-1)
% L2 = uncertainty{1} * A' * temp;

for j = 1:T+1 
    
    energy_RTS(j) = sum( abs( all_states_RTS(:,j) ).^2 );

end

% all_states_RTS(:,1) = all_states_KF(:,1) + (uncertainty{1} * A' ...
%             * (A * uncertainty{1} * A' + Gamma * Q * Gamma')^(-1) ...
%             * (all_states_RTS(:,2) - all_states_KF(:,2)));
% uncertainty_RTS{1} = 

%% 

% comparing the computed solutions

psi = zeros(length(x), length(y));
psi_KF = zeros(length(x), length(y));
psi_RTS = zeros(length(x), length(y));
psi_pred = zeros(length(x), length(y));

which_step = 167;
for j = 1:M*N
        
        % no Stommel solution
        
        psi_KF = psi_KF + all_states_KF(j, which_step) .* exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

        psi_RTS = psi_RTS + all_states_RTS(j, which_step) .* exp(-1i .* beta .* X /vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);
        
        psi = psi + all_states(j, which_step) .* exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

        psi_pred = psi_pred + all_states_pred(j, which_step) .* ...
            exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

end

% adding stommel solution

psi_KF = psi_KF + exp(-X * beta / Ra) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);

psi = psi + exp(-X * beta / Ra) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);

psi_RTS = psi_RTS + exp(-X * beta / Ra) .* sin(pi .* Y) + (X - 1) .* sin(pi .* Y);


%%

% Plotting the modes and their inverses 

% creating my figure 10 

% figure()
% tiledlayout(3,1);
% 
% % non-dimensional periods, frequencies, and the logarithm of the initial
% % condition
% 
% sigma_plot = nexttile;
% plot(2*pi./abs(sigma_nm(:)), '*', 'LineWidth', 1.5);
% ylabel('Non-dim Periods', 'FontSize', 13)
% 
% sigma_inverse = nexttile;
% plot(abs(sigma_nm(:)), '*', 'LineWidth', 1.5);
% ylabel('Radian freq.', 'FontSize', 13)
% 
% log_initial = nexttile;
% semilogy(x0(:), '*', 'LineWidth', 1.5)
% ylabel('Semilog plot of initial cond.', 'FontSize', 13)
% 
% xlabel('1, 2, ..., N * M', 'FontSize', 13)

%%

% just plotting the streamfunction after a few time steps 

figure()

contourf(x,y,real(psi))
one_title = title('$\psi(t, x, y)$', 'FontSize', 19);
one_title.Interpreter = "latex";
xlabel('x', 'FontSize', 13)
ylabel('y', 'FontSize', 13)
colorbar

% add on the locations of observational data
hold on
plot(x_vals, y_vals, 'o', 'linewidth', 4, 'Color', 'r')

% add on where we compute the WBC
hold on
plot(0:.01:0.2, 0.5*ones(size(0:.01:0.2)), 'linewidth', 3);


%%

% plotting the computed solutions 

figure()
tiledlayout(1, 3)

nexttile;
contourf(x,y,real(psi))
title('True Wave')
xlabel('x')
ylabel('y')
colorbar


nexttile;
contourf(x, y, real(psi_KF))
title('KF Wave')
xlabel('x')
ylabel('y')
colorbar

nexttile;
contourf(x, y, real(psi_RTS))
title('RTSWave')
xlabel('x')
ylabel('y')
colorbar


%%

% plotting the "energy" of the system

figure()
tiledlayout(3,1)

t = 0:2000;
nexttile;
plot(t, energy, t, energy_pred, 'linewidth', 1.5)
hold on 
plot(t, energy_KF, '--', 'linewidth', 1.7)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.9) 
xline(data_steps, ':')
legend('True', 'Pred', 'KF', 'RTS', 'FontSize', 12.5)
ylabel('Energy', 'FontSize', 13)

nexttile;
plot(t, energy, t, energy_pred, '--', 'linewidth', 1.5)
hold on 
plot(t, energy_KF, '--', 'linewidth', 1.7)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.7) 
xline(data_steps, ':')
xlim([400, 1350])
legend('True', 'Pred', 'KF', 'RTS', 'FontSize', 12.5)
ylabel('Energy', 'FontSize', 13)

nexttile;
plot(t, energy - energy_KF, t, energy - energy_RTS, '--', 'linewidth', 1.5)
xline(data_steps, ':')
yline(0)
legend('True - KF', 'True - RTS', 'FontSize', 12.5)
ylabel('Energy', 'FontSize', 13)

% nexttile; 
% plot(t, energy_KF, t, energy_RTS, '--', 'linewidth', 1.5)
% xline(data_steps)
% ylabel('Energy', 'FontSize', 13)
% legend('KF', 'RTS')

xlabel('Time step', 'FontSize', 13)

% nexttile;
% plot(energy_KF(1:end))
% title('KF Energy')


%% 

% Here I plot the singular values of E and the diagonal elements of the
% solution resolution matrix 

[U, S, V] = svd(E);

tiledlayout(2,1);

nexttile;
plot(diag(S), '*', 'linewidth', 2)
xlabel('Singular value number')
one = legend('$\lambda_j$', 'FontSize', 15);
one.Interpreter="latex";

nexttile;
plot(diag(V(:, 1:end-6)*V(:, 1:end-6)'), '*', 'linewidth', 2)
xlabel('State vector element')
xticks(1:18)
ylim([0, 1.5])
two = legend('diag$(V V^T)_j$', 'FontSize', 15);
two.Interpreter="latex";

% 
% nexttile;
% plot(diag(V))