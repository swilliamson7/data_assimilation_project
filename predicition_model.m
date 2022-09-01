% This script will create a plot that compares the energy from a
% prediction model to the "true" energy. We're basically running the
% forward model under two cases:
%
% The true model will have the "correct" deterministic component
%                q(t) = 0.1 cos( 2 pi t/(2.5/r) ) + eps(t)
% and no noise in the data points. The prediction model will have the
% "incorrect" deterministic component
%                q_KF(t) = 0.1 * cos(2 pi t / (2.5/r) )
% and there will be noise in the data points.

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

% loading the noise that we'll add to the predicition model
load('noise_0_1_sd.mat')

% create the "truth"
[~, ~, ~, ~, all_states, eps, k_forcing, ...
    kinetic, potential] = forward_func(x0, k, r, dt, M, B, q, Gamma, u, 0.0.*noise);

% now we create the "prediction"
q_KF = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );
[~, ~, ~, ~, all_states_pred, eps_pred, k_forcing_pred, kinetic_pred, ...
                potential_pred] = forward_func(x0, k, r, dt, M, B, q_KF, 0.0*Gamma, 0.0.*u, 0.0.*noise);

t = 0:M-1;            
plot(t, eps, t, eps_pred, '--', 'linewidth', 1.5)
xlabel('Time step')
ylabel('Energy')
legend('"True" energy', 'Prediction energy')