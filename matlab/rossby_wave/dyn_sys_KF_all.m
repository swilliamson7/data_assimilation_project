% This script will compute the coefficients of the modes in our simplified
% dynamical system, not a true solution to the Rossby wave equation. To do 
% this we use the Kalman filter for data assimilation and an RTS smoother
% for smoothing. Wedefine 
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

T = 500;   % how many steps to take 
dt = 29;     % time step
L = 1;       % length of square basin, 1 b/c non-dimensionalized
b = 1.8e-3;
a = 6.3e6;   % radius of Earth, kilometers
theta = 30;  % latitude we use for computation of beta and f 

% Earth's rotation rate
Omega = 7.2921e-5;         % rads / sec 

% determining beta and f at our given latitude
beta = 2 * Omega * cosd(theta) / a;
f = 2 * Omega * sind(theta);

% the choice of n, m determines which of the modes we keep 
n = 3:5;
m = 4:9;
N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  
vec_n = Nn(:);
vec_m = Mm(:);

% initial condition for the system
x0 = (1./(vec_n.^2 + vec_m.^2));

% dispersion relation, given as an N by M matrix, with increases in n
% happening across rows and increases in m happening across columns 
sigma_nm = - beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );

% building the forward operator for stepping the coefficients 

A = exp(-b) .* diag( exp(-1i * sigma_nm(:) * dt ) );    % No time dependence/diss

% random forcing, will change to be consistently the same forcing 
q = @(t) 0.001 * randn(M*N, 1);

% pre-allocating storage for quantities we want to store

% storage for the "true" states and "energy"
state = x0;
all_states = zeros(M*N,T+1);
all_states(:,1) = x0;

energy = zeros(1, T+1);
energy(1) = sum(abs(all_states(:,1)).^2);

% storage for the KF computed states and energy
x0_KF = x0;
all_states_KF = zeros(M*N, T+1);
all_states_KF(:,1) = x0_KF;        % KF has perfect initial condition

energy_KF = zeros(1, T+1);
energy_KF(1) = sum(abs(all_states_KF(:,1)).^2);

% storage for the smoothed states and energy
all_states_RTS = zeros(M*N, T+1);
energy_RTS = zeros(1, T+1);

% noise to add to the data points to make them a little fuzzy

noise = 0.001.*randn(M*N,T+1);

%%

% Running the forward model to create the artificial data
for j = 2:T+1
    
    state = A * state + q(j);

    all_states(:,j) = state;
    
    energy(j) = sum(abs(state).^2);

end

% how to distribute data
E = eye(N * M);

% adding noise here to make the data a little fuzzy
data = E * all_states + noise; 

%%

% Lastly we want to run the Kalman filter and reconstruct the coefficients
% of the normal mode decomposition 

% covariance matrix
uncertainty = cell(T,1);
P0 = zeros(M*N, M*N);
%P0(1,1) = var((1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1));
uncertainty{1} = P0;

% matrix of covariance of noise in the data
R = var(noise(:)) .* eye(N * M);

% covariance of random forcing
Q = var(q(2)) .* eye(N * M);
Gamma = eye(N*M); % we apply forcing to all entries of the state vector

% timesteps where we'll incorporate data
data_steps = 10:10:400;

% running the Kalman filter 

state_old = x0_KF;
P_old = P0;
temp = dt;
for j = 2:T+1

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

% Now that the Kalman filter has been used to the final time and all
% available data has been incorporated, we want to run the RTS smoother
% backwards to improve the prediction 

% uncertainty_RTS = cell(T+1,1);
% 
% for j = T+1:-1:2
% 
%     L = uncertainty{j-1} * A' * uncertainty{j};
%     M = Q * Gamma' * 
% 
% end


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
        
        % with the added stommel solution

        psi_KF = psi_KF + exp(-X * 0.29) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
        
        psi = psi + exp(-X * 0.29) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
        
end

%%

% Plotting the modes and their inverses 

% creating my figure 10 

figure()
tiledlayout(3,1);

% non-dimensional periods

sigma_inverse = nexttile;
plot(1./sigma_nm(:), '*');

sigma_plot = nexttile;
plot(sigma_nm(:), '*');


%%

% plotting the computed solutions 

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