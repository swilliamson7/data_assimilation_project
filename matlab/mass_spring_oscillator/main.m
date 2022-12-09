% here we can play around with the forward function to generate some plots 

x0 = [1;0;0;0;0;0];
k = 30;
r = 0.5;
dt = 0.001;
M = 10000;

%B = 0;
B = x0;
q = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

%Gamma = 0;
Gamma = x0;

load('random_forcing_0_1_sd.mat')

noise = 0.0.*randn(6, M);
[A, Ac, Kc, Rc, all_states, eps, k_forcing, ...
    kinetic, potential] = forward_func(x0, k, r, dt, M, B, q, Gamma, u);

t = 0:M-1;

% Generates figure 2 of the paper

% tiledlayout(4,1)
% 
% x1_x1minx3 = nexttile;
% plot(x1_x1minx3, t, all_states(1,:), t, all_states(1,:) - all_states(3,:), '--', 'linewidth', 1.5)
% legend('x_1', 'x_1 - x_3')
% ylabel('Displacement')
% 
% x4_x5minx6 = nexttile;
% plot(x4_x5minx6, t, all_states(4,:), t, all_states(5,:) - all_states(6,:), '--', 'linewidth', 1.5)
% legend('x_4', 'x_5 - x_6')
% ylabel('Displacement')
% 
% energy = nexttile;
% plot(energy, t, eps, 'linewidth', 1.5)
% ylabel('Energy')
% 
% potential_kinetic_eng = nexttile;
% plot(potential_kinetic_eng, t, potential, t, kinetic, '--', 'linewidth', 1.5)
% ylabel('Energy')
% legend('Kinetic energy', 'Potential energy')
% 
% xlabel('Time step')

% Generates figure 3 of the paper

tiledlayout(3,1)

q_plot = nexttile;
plot(q_plot, t, k_forcing + u)
one = ylabel('$q_1(t)$', 'FontSize', 15);
one.Interpreter = "latex";

x1_x1minx3 = nexttile;
plot(x1_x1minx3, t, all_states(1,:), t, all_states(1,:) - all_states(3, :), '--', 'linewidth', 1.5)
ylabel('Displacement')
two = legend('$x_1(t)$', '$x_1(t) - x_3(t)$', 'FontSize', 14);
two.Interpreter = "latex";

% x4 = nexttile;
% plot(x4, t, all_states(4,:), 'linewidth', 1.5)
% ylabel('Displacement')
% legend('x_4')

tot_energy = nexttile;
plot(tot_energy, t, eps, t, kinetic, '--', t, potential, '-.', 'linewidth', 1.5)
ylabel('Energy')
legend('Total energy', 'Kinetic energy', 'Potential energy', 'FontSize', 14)

xlabel('Time step')
