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

L = 1;  % length of square basin, 1 b/c non-dimensionalized
dx = .01; dy = .01;
x = 0:dx:L; y = 0:dy:L;
[X, Y] = meshgrid(x,y);

T = 2000;    % how many steps to take 
dt = 29;     % time step
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

% inline Stommel solution function 

C1 = exp(1/2 * ( -sqrt(4 * pi^2 * Ra^2 + beta^2)/Ra - beta / Ra ) );
C2 = exp(1/2 * ( sqrt(4 * pi^2 * Ra^2 + beta^2)/Ra - beta / Ra) );

exp_neg = @(x) exp(x/2 .* (-sqrt(4 * pi^2 * Ra^2 + beta^2)/Ra - beta / Ra ));
           
exp_pos = @(x) exp(x/2 .* ( sqrt(4 * pi^2 * Ra^2 + beta^2)/Ra - beta / Ra) );

stommel_soln = @(x, y) 1/(pi^2 * Ra) .* ( ((C2 - 1) / (C2 - C1)) .* exp_neg(x) ...
                + ((1 - C1) / (C2 - C1)) .* exp_pos(x) - 1 ) ...
                .* sin(pi * y);

% initial condition for the system
x0 = (1./(vec_n.^2 + vec_m.^2));

% adding in the one that relates to the Stommel solution
x0(end + 1) = 1;

% dispersion relation, given as an N by M matrix, with increases in m
% happening across rows and increases in n happening across columns 
sigma_nm = beta * L ./ ( 2 * pi * sqrt(Nn.^2 + Mm.^2) );
vec_sigma_nm = sigma_nm(:);

% creating a function that computess a specific sigma_nm value
sigma_func = @(n,m) -beta * L / (2 * pi * sqrt(n^2 + m^2) );

% building the forward operator for stepping the coefficients 
A = exp(-b) .* diag( exp(-1i * sigma_nm(:) * dt) ); 

% adding in the final diagonal element for the Stommel solution 
A(end + 1, end + 1) = 1;

% random forcing
load('new_forcing_1_sd.mat');
q = 0.002 .* forcing;
q(end, :) = 0;

% pre-allocating storage for quantities

% storage for the "true" states and "energy"
all_states = zeros(M*N + 1,T+1);
all_states(:,1) = x0;

energy = zeros(1, T+1);
energy(1) = sum(abs(all_states(1:end-1,1)).^2);

% storage for the KF computed states and energy

x0_KF = 1.5 .* x0;                    % 50percent larger IC for KF  
all_states_KF = zeros(M*N + 1, T+1);
all_states_KF(:,1) = x0_KF;        

% storage for the predicted states and energy
all_states_pred = zeros(M*N + 1, T+1);
all_states_pred(:,1) = x0_KF;

energy_KF = zeros(1, T+1);
energy_KF(1) = sum(abs(all_states_KF(1:end-1,1)).^2);

energy_pred = zeros(1, T+1);
energy_pred(1) = energy_KF(1);

% noise to add to the data points to make them a little fuzzy

load('new_noise_1_sd.mat');
noise = 0.001.*noise;

%%

% Running the forward model to create the artificial data

state = x0;

for j = 2:T+1
    
    % the zero here ensures that we don't add forcing to the steady soln
    state = A * state + q(:,j);

    all_states(:,j) = state;
    
    energy(j) = sum(abs(state(1:end-1)).^2);

end

% now we want to create data from the psi values we just computed. This 
% requires us to build the matrix E differently than in our mass spring 
% oscillator system

x_vals = [1/8:1/8:7/8, 1/8:1/8:7/8];
y_vals = [1/3.*ones(size(1/8:1/8:7/8)), 2/3.*ones(size(1/8:1/8:7/8))];


% x_vals = [.2:.1:.8, .2:.1:.8];
% y_vals = [0.3.*ones(size(.2:.1:.8)), 0.6.*ones(size(.2:.1:.8))];


% if we want to examine the impact of more data use the following for the x
% values and y values 

% x_range = 1/8:1/10:7/8;
% x_vals = [x_range, x_range, x_range, x_range];
% y_vals = [1/5.*ones(size(x_range)), 2/5.*ones(size(x_range)), ...
%           3/5.*ones(size(x_range)), 4/5.*ones(size(x_range))];

E = zeros(length(x_vals), M*N + 1);

for j = 1:length(x_vals)
    for k = 1:M*N

        E(j,k) = exp(-1i * beta * x_vals(j) / vec_sigma_nm(k) ) ...
            * sin(vec_n(k) * pi * x_vals(j)) * sin(vec_m(k) * pi * y_vals(j));

    end

%     E(j,end) = exp(-x_vals(j) * beta / Ra) * sin(pi * y_vals(j)) ...
%                    + (x_vals(j) - 1) * sin(pi * y_vals(j));

      E(j, end) = stommel_soln(x_vals(j), y_vals(j));

end

%E(end, end) = 1;

% adding noise here to make the data a little fuzzy
data = E * all_states + noise(length(x_vals), 1:T+1); 

%%

% Now we want to run the Kalman filter and reconstruct the coefficients
% of the normal mode decomposition 


q_KF = 0.5 .* q;        % 50percent smaller forcing field for KF      
q_KF(end,:) = 0;

% covariance matrix
uncertainty = cell(T+1,1);
P0 = cov( (x0 - x0_KF) * (x0 - x0_KF)' );
uncertainty{1} = P0;

uncertainty_before_data = cell(T+1, 1);

% matrix of covariance of noise in the data
R = var(noise(:)) .* eye(length(x_vals));

% covariance of random forcing
Q = var(q_KF(:)) .* eye(length(x0_KF));
Q(end, end) = 0.0;

% we apply forcing to all entries of the state vector except the steady 
% solution 
Gamma = eye(length(x0_KF)); 
Gamma(end, end) = 0.0;

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
    energy_pred(j) = sum( abs(state_pred(1:end-1)).^2 );

end

states_before_data = zeros(length(x0_KF), T+1);
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
    energy_KF(j) = sum( abs(state_now(1:end-1)).^2 ); 
    
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
all_states_RTS = zeros(length(x0_KF), T+1);
energy_RTS = zeros(1, T+1);

uncertainty_RTS = cell(T+1,1);
Q_RTS = cell(T+1,1);

Q_RTS{end} = Q;

% the initial value for the smoother is just what was computed with the KF
uncertainty_RTS{end} = uncertainty{end};
all_states_RTS(:, end) = all_states_KF(:, end);

forcing_RTS = zeros(length(x0_KF),T+1);
forcing_RTS(:, end) = q_KF(:, end);

L_operators = cell(T, 1);

M_operators = cell(T, 1); 

Q_new = Q;

for j = T:-1:1

    L = uncertainty{j} * A' * uncertainty_before_data{j+1}^(-1);
    M_RTS = Q * Gamma' * uncertainty_before_data{j+1}^(-1);

    all_states_RTS(:, j) = all_states_KF(:, j) + ...
                L * (all_states_RTS(:, j+1) - states_before_data(:, j+1));

    uncertainty_RTS{j} = uncertainty{j} + L * (uncertainty_RTS{j+1} ...
                - uncertainty_before_data{j+1}) * L';

    forcing_RTS(:, j) = q_KF(:, j) + M_RTS * (all_states_RTS(:, j+1) - ...
                                            states_before_data(:, j+1));

    Q_new = Q + M_RTS * (uncertainty_RTS{j+1} - uncertainty_before_data{j+1}) ...
                                * M_RTS';

    Q_RTS{j} = Q_new;
    
    M_operators{j} = M_RTS;
    L_operators{j} = L;

end

for j = 1:T+1 
    
    energy_RTS(j) = sum( abs( all_states_RTS(1:end-1,j) ).^2 );

end

%% 

% creating the psi values for all of the different tools (KF, pred, RTS).

psi = cell(T+1, 1);
psi_KF = cell(T+1, 1);
psi_RTS = cell(T+1, 1);
psi_pred = cell(T+1, 1);

for j = 1:T+1

    psi{j} = zeros(length(x), length(y));
    psi_KF{j} = zeros(length(x), length(y));
    psi_RTS{j} = zeros(length(x), length(y));
    psi_pred{j} = zeros(length(x), length(y));

end

for k = 1:T+1

for j = 1:M*N
        
        % no Stommel solution
        
        psi_KF{k} = psi_KF{k} + all_states_KF(j, k) .* exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

        psi_RTS{k} = psi_RTS{k} + all_states_RTS(j, k) .* exp(-1i .* beta .* X /vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);
        
        psi{k} = psi{k} + all_states(j, k) .* exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

        psi_pred{k} = psi_pred{k} + all_states_pred(j, k) .* ...
            exp(-1i .* beta .* X / vec_sigma_nm(j) ) ...
            .* sin(vec_m(j) .* pi .* Y) .* sin(vec_n(j) .* pi .* X);

end

% adding stommel solution

% psi_KF{k} = psi_KF{k} + all_states_KF(end,k) .* (exp(-X * beta / Ra) .* ...
%         sin(pi .* Y) + (X - 1).* sin(pi .* Y));

psi_KF{k} = psi_KF{k} + all_states_KF(end, k) .* stommel_soln(X, Y);

% psi{k} = psi{k} + all_states(end, k) .* (exp(-X * beta / Ra) .* ...
%         sin(pi .* Y) + (X - 1).* sin(pi .* Y));
psi{k} = psi{k} + all_states(end, k) .* stommel_soln(X, Y);

% psi_RTS{k} = psi_RTS{k} + all_states_RTS(end,k) .* (exp(-X * beta / Ra) ...
%         .* sin(pi .* Y) + (X - 1) .* sin(pi .* Y));
psi_RTS{k} = psi_RTS{k} + all_states_RTS(end,k) .* stommel_soln(X, Y);

% psi_pred{k} = psi_pred{k} + all_states_pred(end,k) .* (exp(-X * beta / Ra) ...
%         .* sin(pi .* Y) + (X - 1) .* sin(pi .* Y));

psi_pred{k} = psi_pred{k} + all_states_pred(end,k) .* stommel_soln(X, Y);

end

% psi_stommel = exp(-X * beta / Ra) .* sin(pi .* Y) + (X - 1) .* sin(pi .* Y);
psi_stommel = stommel_soln(X, Y);

%% 

% here we compute the western boundary current from the computed psi values

wbc_true = zeros(T+1, 1);
wbc_pred = zeros(T+1, 1);
wbc_KF = zeros(T+1, 1);
wbc_RTS = zeros(T+1, 1);

wbc_stommel = zeros(T+1,1);

for j = 1:T+1

    wbc_true(j) = real(psi{j}(1, 50)) - real(psi{j}(20, 50));

    wbc_pred(j) = real(psi_pred{j}(1, 50)) - real(psi_pred{j}(20, 50));

    wbc_KF(j) = real(psi_KF{j}(1, 50)) - real(psi_KF{j}(20, 50));

    wbc_RTS(j) = real(psi_RTS{j}(1, 50)) - real(psi_RTS{j}(20, 50));

    wbc_stommel(j) = real(psi_stommel(1, 50)) - real(psi_stommel(20,50));

end
%%

% computing and storing the uncertainties wanted for plots

steady_uncertainty_KF = zeros(1, T+1);
steady_uncertainty_RTS = zeros(1, T+1); 

L_matrix_diag_real = zeros(1, T);
L_matrix_diag_imag = zeros(1, T);
L_norms = zeros(1, T);

M_norms = zeros(1, T); 
M_multiplied_norms = zeros(1,T);

forcing_RTS_uncertainty = zeros(1,T+1);

for j = 1:T+1

    steady_uncertainty_KF(j) = sqrt(real(uncertainty{j}(end, end)));
    steady_uncertainty_RTS(j) = sqrt(real(uncertainty_RTS{j}(end, end)));

    forcing_RTS_uncertainty(j) = Q_RTS{j}(2,2)^(1/2);

end

for j = 1:T

    L_matrix_diag_real(j) = real(L_operators{j}(2,2));
    L_matrix_diag_imag(j) = imag(L_operators{j}(2,2)); 
    L_norms(j) = norm(L_operators{j});

    M_norms(j) = norm(M_operators{j});

end

M_multiplied_norms(T) = norm(M_operators{j});
temp = M_operators{T};
for j = T-1:-1:1
    
    temp = M_operators{j} * temp;
    M_multiplied_norms(j) = norm(temp);

end

%% computing the vector that we can use to find the uncertainty of the WBC 

F = zeros(N * M + 1, 1);

for j = 1:M*N

F(j) = exp(-1i .* beta .* 0 /vec_sigma_nm(j) ) ...
            * sin(vec_m(j) .* pi .* 0.5) .* sin(vec_n(j) .* pi .* 0) ...
       - exp(-1i .* beta .* 0.2 /vec_sigma_nm(j) ) ... 
            * sin(vec_m(j) .* pi .* 0.5) .* sin(vec_n(j) .* pi .* 0.2);

end

% F(end) =  exp(-0 * beta / Ra) * sin(pi .* 0.5) + (0 - 1) .* sin(pi .* 0.5)...
%         - exp(-0.2 * beta / Ra) * sin(pi .* 0.5) + (0.2 - 1) .* sin(pi .* 0.5);

F(end) = stommel_soln(0, 0.5) - stommel_soln(0.2, 0.5);

%% using F to find the uncertainty values in the WBC estimates 

wbc_uncertainty_KF = zeros(T+1, 1);
wbc_uncertainty_RTS = zeros(T+1, 1);

for j = 1:T+1 
    
    wbc_uncertainty_KF(j) = sqrt(F' * uncertainty{j} * F);

    wbc_uncertainty_RTS(j) = sqrt(F' * uncertainty_RTS{j} * F);

end


%% Generates the plots for the rossby wave


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything past this point is a plot of something 
t = 0:2000;

tiledlayout(3,1);

nexttile;
plot(t, wbc_KF)
hold on

% fill([t, fliplr(t)], [wbc_KF' + wbc_uncertainty_KF', fliplr(wbc_KF' - wbc_uncertainty_KF')], ...
%                         'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.3)
% hold on
% nexttile;
% shadedErrorBar(t, wbc_KF, wbc_uncertainty_KF, 'lineprops','-b','patchSaturation',0.000001);
% hold on 
plot(t, wbc_stommel, '--', 'linewidth', 1.5)
hold on
xline(data_steps, ':')
%xlabel('Time step', 'FontSize', 14)
legend('KF', 'Steady solution', 'FontSize', 15)



nexttile;
plot(t, wbc_RTS)
hold on
plot(t, wbc_stommel, '--', 'linewidth', 1.5)
xline(data_steps, ':')
%xlabel('Time step', 'FontSize', 14)
legend('RTS', 'Steady solution', 'FontSize', 15)

nexttile; 

plot(t, wbc_uncertainty_KF, t, wbc_uncertainty_RTS, '--', 'linewidth', 1.5)
xline(data_steps, ':')
legend('KF Uncertainty', 'RTS Uncertainty', 'fontsize', 15)
xlabel('Time step', 'FontSize', 14)


%%

% Plotting the modes and their inverses 

figure()
tiledlayout(2,1);

% non-dimensional periods, frequencies, and the logarithm of the initial
% condition

sigma_plot = nexttile;
plot(2*pi./abs(sigma_nm(:)), '*', 'LineWidth', 1.5);
hold on
yline(dt, '--', 'LineWidth', 2)
ylabel('Non-dim Periods', 'FontSize', 13)

sigma_inverse = nexttile;
plot(abs(sigma_nm(:)), '*', 'LineWidth', 1.5);
ylabel('Radian freq.', 'FontSize', 13)

% log_initial = nexttile;
% semilogy(x0(:), '*', 'LineWidth', 1.5)
% ylabel('Semilog plot of initial cond.', 'FontSize', 13)

xlabel('1, 2, ..., N * M', 'FontSize', 13)

%%

% just plotting the streamfunction after a few time steps 

figure()

contourf(x,y,real(psi{167}))
one_title = title('$\psi(167 \Delta t, x, y)$', 'FontSize', 19);
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
contourf(x,y,real(psi{167}))
title('True Wave')
xlabel('x')
ylabel('y')
colorbar

nexttile;
contourf(x, y, real(psi_KF{167}))
title('KF Wave')
xlabel('x')
ylabel('y')
colorbar

nexttile;
contourf(x, y, real(psi_RTS{167}))
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
%plot(t, energy, 'linewidth', 1.5)
hold on 
plot(t, energy_KF, '--', 'linewidth', 1.9)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.5) 
xline(data_steps, ':')
legend('True', 'Pred', 'KF', 'RTS', 'FontSize', 12.5)
%legend('True', 'KF', 'RTS', 'FontSize', 12.5)
ylabel('Energy', 'FontSize', 13)

nexttile;
plot(t, energy, t, energy_pred, '--', 'linewidth', 1.5)
%plot(t, energy, 'linewidth', 1.5)
hold on 
plot(t, energy_KF, '--', 'linewidth', 1.9)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.5) 
xline(data_steps, ':')
xlim([400, 1350])
legend('True', 'Pred', 'KF', 'RTS', 'FontSize', 12.5)
%legend('True', 'KF', 'RTS', 'FontSize', 12.5)
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


%%
% 
% % plotting the average of psi over all timesteps 
% 
% average_psi = zeros(101, 101);
% 
% for j = 2:2001 
%     
%     average_psi = average_psi + psi_KF{j} / 2001;
% 
% end
% 
% tiledlayout(1,2);
% 
% nexttile;
% contourf(x, y, real(average_psi))
% xlabel('x', 'FontSize', 13)
% ylabel('y', 'FontSize', 13)
% title('Averaged KF')
% colorbar
% 
% nexttile;
% contourf(x, y, psi_stommel)
% xlabel('x', 'FontSize', 13)
% ylabel('y', 'FontSize', 13)
% title('Stommel solution')
% colorbar


%%

tiledlayout(2,1);

nexttile;
plot(t(2:end), L_norms, 'LineWidth', 1.5)
xline(data_steps, ':')
xlabel('Time step', 'FontSize', 14)
two = legend('$||\mathbf{L}(t)||$', 'FontSize', 14);
two.Interpreter = "latex";

nexttile; 
plot(t(2:end), L_matrix_diag_real, t(1:end-1), L_matrix_diag_imag, '--', ...
            'LineWidth', 1.5)
xline(data_steps, ':')
one = legend('Re$(L_{22}(t))$', 'Im$(L_{22}(t))$', 'FontSize', 14);
one.Interpreter = "latex";

xlabel('Time step', 'FontSize', 14)

%%
tiledlayout(1,1);

nexttile;
plot(t(1:end-1), M_norms, 'linewidth', 1.5)
xline(data_steps, ':')
xlabel('Time step', 'FontSize', 14)
one = legend('$||\mathbf{M}(t)||$', 'FontSize', 17);
one.Interpreter = "latex";
% 
% nexttile;
% plot(-20:1:0, M_multiplied_norms(end:-1:end-20))

% nexttile;
% plot(t, steady_uncertainty_KF, t, steady_uncertainty_RTS, '--', 'linewidth', 1.5)
% xline(data_steps, ':')
% legend('KF Uncertainty', 'RTS Uncertainty', 'FontSize', 14)
% xlabel('Time step')


%%

tiledlayout(2,1);

nexttile;
plot(t, forcing_RTS(2, :))
forcing = legend('$u_2(t)$', 'FontSize', 17);
forcing.Interpreter = "latex";

nexttile;
plot(t, forcing_RTS_uncertainty, 'linewidth', 1.5)
xline(data_steps, ':')
two = legend('$(\mathbf{Q}_{2, 2}(t, +))^{1/2}$', 'FontSize', 17);
two.Interpreter = "latex";
xlabel('Timestep', 'FontSize', 15)

%% 
% 
% tiledlayout(2,1);
% 
% nexttile;
% plot(t, real(all_states_KF(end, :)), 'linewidth', 1.5)
% hold on
% std_dev_plot = errorbar(t, real(all_states_KF(end, :)), steady_uncertainty_KF);
% alpha = 0.1;   
% % Set transparency (undocumented)
% set([std_dev_plot.Bar, std_dev_plot.Line], 'ColorType', 'truecoloralpha', 'ColorData', [std_dev_plot.Line.ColorData(1:3); 100*alpha])
% 
% nexttile;
% plot(t, real(all_states_RTS(end, :)), 'linewidth', 1.5)
% hold on
% std_dev_plot = errorbar(t, real(all_states_RTS(end, :)), steady_uncertainty_RTS);
% alpha = 0.1;   
% % Set transparency (undocumented)
% set([std_dev_plot.Bar, std_dev_plot.Line], 'ColorType', 'truecoloralpha', 'ColorData', [std_dev_plot.Line.ColorData(1:3); 100*alpha])
% hold on