% We now want to apply the KF to a forced Rossby wave system. This script
% will use the function rossby_forward_forced which runs the system with
% random forcing at every step. We'll assume that the KF does not know
% about the forcing, but has perfect knowledge of A, and near-perfect
% knowledge of the initial condition

% the solution used here is for periodic x-boundary and zero y-boundary, so
% not the solution that we wanted to discuss in the manuscript 

clc, clear

% deciding which modes to keep
n = 1:5;
m = 1:6;
N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  
vec_n = Nn(:);
vec_m = Mm(:);

beta = 1.7;
L = 1;

dt = 0.01;

% a function for the dispersion relation, depends on n and m, needed to
% calculate energy later

sigma_func = @(n, m) -beta * (L / 2) / ( sqrt(n^2 * pi^2 + m^2 * pi^2) );

% the number of time steps to take
T = 10000;

% random forcing
q_nm = @(t) 0.1 * randn(N*M, 1);

% constant initial condition, chosen so that the sum over normal modes
% converges
x0 = 1./(vec_n.^2 + vec_m.^2);

% noise for the data
noise = 0.0001.*randn(N*M, T+1);

% creating the data for the system
[A, sigma_nm, all_states, energy] = rw_forward_func_forced(x0, dt, T, n, m, 1, 1, q_nm);

% now we want to run the KF and see if we can reconstruct the coefficients

data = all_states + noise;

% where to apply data points

dataset = 1000:75:3000;
%dataset = 1:T+1; 

% initial condition (and thus uncertainty) seen by KF
%x0_KF = x0 + (1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1);
x0_KF = x0;
P0 = zeros(M*N, M*N);
P0(1,1) = var((1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1));

% how to distribute data
E = eye(N*M);

% matrix of covariance of noise in the data
R = var(noise(:)) .* eye(N*M);

% covariance of noise in the forcing
Q = var(q_nm(2)) .* eye(N * M);
Gamma = eye(N*M); % we apply forcing to all entries of the state vector

% storing the states found by the KF
all_states_KF = zeros(N*M, T+1);
all_states_KF(:, 1) = x0_KF;

% place to store covariance matrix P
uncertainty = cell(N*M,T+1);
uncertainty{1} = P0;

% place to store energy
energy_KF = zeros(1, T+1);

state_old = x0_KF;
energy_KF(1) = sum(abs(state_old).^2);
P_old = P0;
temp = dt;
for j = 2:T+1
    
    % step forward
    %temp_state = A * state_old;
    temp_state = state_old + A * state_old;
    temp_P = A * P_old * A' + Gamma * Q * Gamma';
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == dataset)
        
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
        
%         psi = psi + mode_coeff(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
%             .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y) ...
%             + exp(-X) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
        

end

%%

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
xline(dataset*dt)

% nexttile;
% plot(energy_KF(1:end))
% title('KF Energy')



