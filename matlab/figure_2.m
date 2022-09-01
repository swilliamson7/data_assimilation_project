% This script will generate the updated Figure 2. 
% In this I'm just running the forward model and plotting the same
% quantities that Carl has in his paper. There is no random forcing or
% deterministic forcing, and no KF. 

clc, clear 

dt = 0.001;
M = 10000;

% initial condition
x0 = [1;0;0;0;0;0];

% spring coefficient and dissipation constant
k = 30;
r = 0.5;

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

state_old = x0;
for j = 2:M
    
    % step forward
    state_now = A * state_old;

    % calculate and store system energy
    kinetic(j) = 0.5 * ( state_now(4:6)' * state_now(4:6) );
    potential(j) = 0.5 * ( - state_now(1:3)' * Kc * state_now(1:3) );
    eps(j) = kinetic(j) + potential(j);

    % store new state
    all_states(:,j) = state_now;

    % cycle states
    state_old = state_now;

end 

t = 1:10000;

tiledlayout(4,1)

x1_x1minx3 = nexttile;
plot(x1_x1minx3, t, all_states(1,:), t, all_states(1,:) - all_states(3,:), '--', 'linewidth', 1.5)
legend('x_1', 'x_1 - x_3', 'Location', 'best')
ylabel('Displacement')

x4_x5minx6 = nexttile;
plot(x4_x5minx6, t, all_states(4,:), t, all_states(5,:) - all_states(6,:), '--', 'linewidth', 1.5)
legend('x_4', 'x_5 - x_6')
ylabel('Displacement')

energy = nexttile;
plot(energy, t, eps, 'linewidth', 1.5)
ylabel('Energy')

potential_kinetic_eng = nexttile;
plot(potential_kinetic_eng, t, kinetic, t, potential, '--', 'linewidth', 1.5)
ylabel('Energy')
legend('Kinetic energy', 'Potential energy')

xlabel('Time step')
