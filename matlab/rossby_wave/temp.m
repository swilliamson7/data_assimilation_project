clc, clear 

n = 3:5;
m = 4:9;
N = length(n);
M = length(m);
T = 2000;

noise = 0.001 .* randn(M*N, T+1);
save('noise_0_001_sd.mat', 'noise')

noise = randn(M*N, T+1);
save('noise_1_sd.mat', 'noise')

forcing = randn(M*N, T+1);
save('forcing_1_sd', 'forcing')
