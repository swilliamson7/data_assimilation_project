% This script will compute the coefficients of the modes in our simplified
% dynamical system, not a true solution to the Rossby wave equation. To do 
% this we use the Kalman filter for data assimilation. We
% define 
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

%%

clc, clear 

T = 10000;   % how many steps to take 
dt = 29;
L = 1;
beta = 1.7;

% the choice of n, m determines which of the modes we keep 
n = 1:5;
m = 1:6;
N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  
vec_n = Nn(:);
vec_m = Mm(:);

% initial condition for the system
x0 = 1./(vec_n.^2 + vec_m.^2);

% dispersion relation, given as an N by M matrix, with increases in n
% happening across rows and increases in m happening across columns 
sigma_nm = - beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );

% building the forward operator for stepping the coefficients 
A = diag( exp(-1i * sigma_nm(:) * dt ) );

% random forcing, will change to be consistently the same forcing 
q = @(t) 0.01 * randn(M*N, 1);

% pre-allocating storage for quantities we want to store, first for the 
% "true" state and then for the Kalman filter results
state = x0;
all_states = zeros(M*N,T+1);
all_states(:,1) = x0;

energy = zeros(1, T+1);
energy(1) = sum(abs(all_states(:,1)).^2);

x0_KF = x0;
all_states_KF = zeros(M*N, T+1);
all_states_KF(:,1) = x0_KF;        % KF has perfect initial condition

energy_KF = zeros(1, T+1);
energy_KF(1) = sum(abs(all_states_KF(:,1)).^2);

% noise to add to the data points to make them a little fuzzy

noise = 0.001.*randn(M*N,T+1);

%%

% building operators used in the Kalman filter

% covariance matrix
uncertainty = cell(M,1);
P0 = zeros(M*N, M*N);
%P0(1,1) = var((1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1));
uncertainty{1} = P0;

% how to distribute data
E = eye(N * M);

% matrix of covariance of noise in the data
R = var(noise(:)) .* eye(N * M);

% covariance of noise in the forcing
Q = var(q(2)) .* eye(N * M);
Gamma = eye(N*M); % we apply forcing to all entries of the state vector

%%

% Now we want to run the forward model to create the artificial data

for j = 2:T+1

    state = A * state + q(j);
    all_states(:,j) = state;
    
    energy(j) = sum(abs(state).^2);

end

% adding noise here to make the data a little fuzzy
data = E * all_states + noise; 

%%

% Lastly we want to run the Kalman filter and reconstruct the coefficients
% of the normal mode decomposition 

% timesteps where we'll incorporate data
data_steps = 1000:75:3000;

% running the Kalman filter 

state_old = x0_KF;
P_old = P0;
temp = dt;
for j = 2:T+1
    
    % step forward
    temp_state = A * state_old;
    temp_P = A * P_old * A' + Gamma * Q * Gamma';
    
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

% comparing the computed solutions

L = 1;
dx = .01; dy = .01;
x = 0:dx:L; y = 0:dy:L;
[X, Y] = meshgrid(x,y);

vec_sigma_nm = sigma_nm(:);

psi = zeros(length(x), length(y));
psi_KF = zeros(length(x), length(y));

for j = 1:M*N
        
        % no Stommel solution
        
        psi_KF = psi_KF + all_states_KF(j, T) .* exp(-1i .* pi .* X .* vec_n(j) ) ...
            .* sin(vec_m(j) .* pi .* Y);
        
        psi = psi + all_states(j, T) .* exp(-1i .* pi .* X .* vec_n(j) ) ...
            .* sin(vec_m(j) .* pi .* Y);
        
%         % with the added stommel solution
% 
%         psi_KF = psi_KF + exp(-X * 0.29) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
%         
%         psi = psi + exp(-X * 0.29) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
%         

end

%%

% here we can compare the computed solutions 

figure()
tiledlayout(1,2)

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

figure()

nexttile;
plot(0:dt:T*dt, energy, 0:dt:T*dt, energy_KF, '--', 'linewidth', 1.5)
title('True Energy')
xline(data_steps*dt)

% nexttile;
% plot(energy_KF(1:end))
% title('KF Energy')
