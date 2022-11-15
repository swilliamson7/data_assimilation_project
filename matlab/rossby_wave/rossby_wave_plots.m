% Generates the plots for the rossby wave


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything past this point is a plot of something 
t = 0:2000;

tiledlayout(2,1);

% nexttile;
% plot(t, wbc_true)
% hold onhm  
% plot(t, wbc_stommel, '--', 'linewidth', 4)

nexttile;
plot(t, wbc_KF)
hold on
plot(t, wbc_stommel, '--', 'linewidth', 1.5)
ylim([0,0.6])
xline(data_steps, ':')
xlabel('Time step', 'FontSize', 14)
legend('KF', 'Steady solution', 'FontSize', 15)

nexttile;
plot(t, wbc_RTS)
hold on
plot(t, wbc_stommel, '--', 'linewidth', 1.5)
ylim([0,0.6])
xline(data_steps, ':')
xlabel('Time step', 'FontSize', 14)
legend('RTS', 'Steady solution', 'FontSize', 15)


%%

% Plotting the modes and their inverses 

% creating my figure 10 

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
plot(t, energy_KF, '--', 'linewidth', 1.7)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.9) 
xline(data_steps, ':')
legend('True', 'Pred', 'KF', 'RTS', 'FontSize', 12.5)
%legend('True', 'KF', 'RTS', 'FontSize', 12.5)
ylabel('Energy', 'FontSize', 13)

nexttile;
plot(t, energy, t, energy_pred, '--', 'linewidth', 1.5)
%plot(t, energy, 'linewidth', 1.5)
hold on 
plot(t, energy_KF, '--', 'linewidth', 1.7)
hold on
plot(t, energy_RTS, ' -.', 'linewidth', 1.7) 
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

% random last minute plots 

tiledlayout(1,1);

nexttile;
plot(t, steady_uncertainty_KF, t, steady_uncertainty_RTS, '--', 'linewidth', 1.5)
xline(data_steps, ':')
legend('KF Uncertainty', 'RTS Uncertainty', 'FontSize', 14)
xlabel('Time step')


%%

tiledlayout(2,1);

nexttile;
plot(t(2:end), L_norms, 'LineWidth', 1.5)
xline(data_steps, ':')
xlabel('Time step', 'FontSize', 14)

nexttile; 
plot(t(2:end), L_matrix_diag_real, t(1:end-1), L_matrix_diag_imag, '--', ...
            'LineWidth', 1.5)
xline(data_steps, ':')
one = legend('Re$(L_{22}(t))$', 'Im$(L_{22}(t))$', 'FontSize', 14);
one.Interpreter = "latex";

xlabel('Time step', 'FontSize', 14)

%%

plot(t(2:end), M_norms)
xline(data_steps, ':')


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