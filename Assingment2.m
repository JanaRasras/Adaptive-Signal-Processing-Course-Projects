%%
% Jana Rusrus
% Assignment 2: FIR causal Wiener solution
%% 
clear all
close all
clc
%
% pkg load signal

%% Question 2

%------------------ setup A ------------------%
% W opt 
N = 10;
n = -N:N-1;
u =@(l) l>=0;

sysA.p2.phi_xx =@(l) 25.25   * (0.98.^l .* u(l) + 0.98.^(-l).* u(-l-1));
sysA.p2.phi_dx =@(l) 824.834 * 0.98.^l - 799.584 * 0.95.^l .* u(l) + ...
             366.125 * 0.98.^(-l).* u(-l-1) + 340.875 * 0.95 .^l .* u(l);

sysA.p2.segma_d = 7249.25;

sysA.p2.P = sysA.p2.phi_dx(0 : N-1)';
sysA.p2.R = toeplitz(sysA.p2.phi_xx(0: N-1)');

sysA.p2.W_opt = sysA.p2.R\sysA.p2.P;
sysA.p2.norm_MMSE = 1 - ( conj(sysA.p2.P')/sysA.p2.R*sysA.p2.P ) / sysA.p2.segma_d

figure()
stem(sysA.p2.phi_dx(0 : N-1)')
input('?')
%subplot(211); stem(n, 0.95.^n .* u(n), 'filled', 'MarkerSize',3)
%title('h[n]'); axis([-1, N, 0, 1.1])

%subplot(212); 
stem(n,  sysA.p2.W_opt', 'filled', 'MarkerSize',3)
title('w[n]Theoritcal_ System A'); axis([-1, N, 0, 10])

%------------------ setup B ------------------%
n = 0:N-1;

sysB.p2.phi_xx =@(l) 25.25 * (0.98.^l .* u(l) + 0.98.^(-l).* u(-l-1));
sysB.p2.phi_dx =@(l) 25.25 * (0.98.^(l+1) .* u(l+1) + 0.98.^(-l-1).* u(-l-2));

sysB.p2.P = sysB.p2.phi_dx(0 : N-1)';
sysB.p2.R = toeplitz(sysB.p2.phi_xx(0: N-1)');
sysB.p2.segma_d = 25.25;

sysB.p2.W_opt = sysB.p2.R\sysB.p2.P;
sysB.p2.norm_MMSE = 1 - ( conj(sysB.p2.P')/sysB.p2.R*sysB.p2.P ) / sysB.p2.segma_d

figure()
%subplot(211); stem(n, 0.95.^n .* u(n), 'filled', 'MarkerSize',3)
%title('h[n]'); axis([-1, N, 0, 1.1])

%subplot(212); 
stem(n, sysB.p2.W_opt', 'filled', 'MarkerSize',3)
title('w[n]Theoritcal_ System B'); axis([-1, N, 0, 2])

%% Question 3
N = 10;
L = 10000000;
u = randn(L+200,1); % white Gaussian, variance of one
x = filter(1,[1 -0.98],u); % t(n) filter
x = x(101:end); % discard first 100 samples to remove transients in simulated signals
v = randn(L+100,1); % white Gaussian, variance of one

% setup A
dA = filter(1,[1 -0.95],x) + v; % h(n) filter plus additive noise
dA = dA(101:end); % discard first 100 samples to remove transients in simulated signals
xA = x(101:end); % discard first 100 samples to remove transients in simulated signals

% setup B
dB = x(101:end); % discard first 100 samples to remove transients in simulated signals
xB = x(100:end-1); % x(n-1), delayed by one sample

%------------------ setup A ------------------%

[phi_est_xx,lagsx] = xcorr(xA, xA,'unbiased');
[phi_est_dx,lagsd] = xcorr(dA, xA,'unbiased');
[phi_est_dd,lagsdd] = xcorr(dA, dA,'unbiased');

% W opt 
n1 = find(lagsx == 0);
n2 = find(lagsd == 0);
n3 = find(lagsdd == 0);

P1_est = phi_est_dx(n2:n2+N-1);
R1_est = toeplitz(phi_est_xx(n1:n1+N-1));

W1_opt_est = R1_est\P1_est;

% MMSE
norm_MMSE1_est = 1 - ( conj(P1_est')/R1_est * P1_est ) / phi_est_dd(n3)

% plotting W_opt
n = 0:N-1;
figure(); 
%subplot(211); stem(n,  W_opt', 'filled', 'MarkerSize',3)
%title('w[n]Theoretical_ System A'); axis([-1, N, 0, 10])

%subplot(212);
stem(n,  W1_opt_est', 'filled', 'MarkerSize',3)
title('w[n]_ est _ System A'); axis([-1, N, 0, 10])

%------------------ setup B ------------------%
[phi2_est_xx,lagsx2] = xcorr(xB, xB,'unbiased');
[phi2_est_dx,lagsd2] = xcorr(dB, xB,'unbiased');
[phi2_est_dd,lagsdd2] = xcorr(dB, dB,'unbiased');

% W opt 
n4 = find(lagsx2 == 0);
n5 = find(lagsd2 == 0);
n6 = find(lagsdd2 == 0);

P2_est = phi2_est_dx(n5:n5+N-1);
R2_est = toeplitz(phi2_est_xx(n4:n4+N-1));

W2_opt_est = R2_est\P2_est;

%------------------ MMSE ------------------%
norm_MMSE2_est = 1 - ( conj(P2_est')/R2_est * P2_est ) / phi2_est_dd(n6)

n = 0:N-1;
figure(); stem(n,  W2_opt_est', 'filled', 'MarkerSize',3)
title('w[n]_ est_ system B'); axis([-1, N, 0, 2])

%% PLOTTING
%------------------ setup A ------------------%
figure()
subplot(211); 
stem(n,  sysA.p2.W_opt', 'filled', 'MarkerSize',3)
title('w[n]Theoritcal_ System A'); axis([-1, N, 0, 10])

subplot(212);
stem(n,  W1_opt_est', 'filled', 'MarkerSize',3)
title('w[n]_ est _ System A'); axis([-1, N, 0, 10])


%------------------ setup B ------------------%
figure()
subplot(211); 
stem(n, sysB.p2.W_opt', 'filled', 'MarkerSize',3)
title('w[n]Theoritcal_ System B'); axis([-1, N, 0, 2])

subplot(212); 
stem(n,  W2_opt_est', 'filled', 'MarkerSize',3)
title('w[n]_ est_ system B'); axis([-1, N, 0, 2])

