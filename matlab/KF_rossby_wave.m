clc, clear

% this script aims to apply the KF to the rossby wave system to reconstruct
% the coeffecients of the normal modes. This **does not** include forcing
% and will be changed

% deciding which modes to keep
n = 1:5;
m = 1:6;
N = length(n); M = length(m);
[Nn,Mm] = meshgrid(n',m');  
vec_n = Nn(:);
vec_m = Mm(:);

beta = 1;
L = 1;

dt = 120;

% a function for the dispertion relation, depends on n and m, needed to
% calculate energy later

sigma_func = @(n, m) -beta * (L / 2) / ( sqrt(n^2 * pi^2 + m^2 * pi^2) );

% the number of time steps to take
T = 20;
t = 0:dt:20;

% constant initial condition, chosen so that the sum over normal modes
% converges
x0 = 1./(vec_n.^2 + vec_m.^2);

% noise for the data
noise = 0.01.*randn(1, T+1);

forcing = 0.5.*randn(1,T+1);
Q = zeros(N*M, N*M);

% creating the data for the system
[A, sigma_nm, all_states, energy] = rossby_wave_forward_func(x0, dt, T, n, m, 1, 1);

% now we want to run the KF and see if we can reconstruct the coefficients

data = all_states + noise;

% where to apply data points

%dataset_1 = [3, 4, 5]; dataset_2 = [10];
dataset_1 = 1:T+1; dataset_2 = dataset_1;

% initial condition (and thus uncertainty) seen by KF
x0_KF = x0 + (1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1);
P0 = zeros(M*N, M*N);
P0(1,1) = var((1./(vec_m.^2 + vec_n.^2)).*0.1.*randn(1));

% how to distribute data
E = eye(N*M);

% matrix of covariance of noise in the data
R = var(noise) .* eye(N*M);

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
    temp_state = A * state_old;
    temp_P = A * P_old * A';
    
    % calculate Kalman matrix
    K = temp_P * E' * ( E * temp_P * E' + R )^(-1);
    
    if sum(j == dataset_1) || sum(j == dataset_2)
        
        state_now = temp_state + K * (data(:, j) - E * temp_state );
        P_now = temp_P - K * E * temp_P;
    else
        state_now = temp_state;
        P_now = temp_P;
        
    end
    
    % computing "energy," not quite sure I agree with this definition of
    % energy
    energy_KF(j) = sum(abs(state_now).^2); 
    
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
        
        psi_KF = psi_KF + all_states_KF(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
            .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y);
        
        psi = psi + all_states(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
            .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y);
        
        % with the added stommel solution
        
%         psi = psi + mode_coeff(j, T) .* exp(-1i .* X ./ vec_sigma_nm(j) ) ...
%             .* sin(vec_n(j) .* pi .* X) .* sin(vec_m(j) .* pi .* Y) ...
%             + exp(-X) .* sin(pi .* Y) + (X - 1).* sin(pi .* Y);
        

end

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
tiledlayout(1,2)

nexttile;
plot(energy)
title('True Energy')
ylim([0.4, 0.41])

nexttile;
plot(energy_KF)
title('KF Energy')



