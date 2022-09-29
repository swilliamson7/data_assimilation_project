clc, clear

x0 = [1;0;0;0;0;0];
k = 30;
r = 0.5;
dt = 0.001;
M = 10000;

% known forcing, only applied to position of mass one
B = x0;
q = @(t) 0.1 * cos( 2 * pi * t / (2.5 / r) );

Gamma = B;
load('random_forcing_0_1_sd.mat')

[Ac, Kc, Rc] = build_matrices(k,r);

% build discrete forward operator 

A = eye(6) + dt.*Ac;

% place to store state vector
all_states = zeros(6,M);
all_states(:,1) = x0;

% place to store energy of the system
eps = zeros(1,M);
kinetic = zeros(1,M);
potential = zeros(1,M);
kinetic(1) = 0.5 * ( x0(4:6)' * x0(4:6) );
potential(1) = 0.5 * (- x0(1:3)' * Kc * x0(1:3) );
eps(1) = kinetic(1) + potential(1);

% place to store forcing in the system 
forcing = zeros(1,M); 
k_forcing = zeros(1,M);

state_old = x0;
stuff = 0;
for j = 2:M
    
    k_forcing(j) = q(stuff);
    forcing(j) = k_forcing(j) + u(j);
    
    % step forward
    state_now = A * state_old + B * k_forcing(j) + B * u(j);

    % calculate and store system energy
    kinetic(j) = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential(j) = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps(j) = kinetic(j) + potential(j);

    % store new state
    all_states(:,j) = state_now;

    % cycle states
    state_old = state_now;
    
    stuff = stuff + dt;
    
end 

tiledlayout(4,1);

t = 0:M-1;


nexttile;
plot(t, forcing)
ylabel('q(t)')

nexttile;
plot(t, all_states(1,:), t, all_states(1,:) - all_states(3,:), '--', 'linewidth', 1.5)
ylabel('Displacement')
one = legend('$x_1(t)$', '$x_1(t) - x_3(t)$', 'FontSize', 13);
one.Interpreter = "latex";

nexttile;
plot(t, eps, 'linewidth', 1.5)
two = legend('$\mathcal{E}(t)$', 'FontSize', 13);
two.Interpreter = "latex";
ylabel('Energy')

nexttile;
plot(t, kinetic, t, potential, '--', 'linewidth', 1.5)
ylabel('Energy')
legend('Kinetic energy', 'Potential energy', 'FontSize', 13)
